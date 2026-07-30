// dh_outlets -- the shared outlet re-derivation for the tiled/distributed DepressionHierarchy builds
// (ENH-5), used by the in-process stitch (tools/dephier_stitch.cpp) and BOTH scans in the distributed
// harness (tools/dephier_mpi.cpp: the full-grid oracle and the per-rank distributed path).
//
// WHY a re-derivation at all: the tiled build discards FloodAndAssignDepressions's per-tile `tile.outlets` and rebuilds
// the outlet set from the RESOLVED (post-conduit) label grid, because BOUNDARY endpoints are only
// meaningful after conduit resolution and each rank must work from its own resolved labels (see the
// ENH-2 rationale in dephier_stitch.cpp). That scan was independently TRIPLICATED, and the same
// NoData-as-ocean bug (a basin->NoData-ocean adjacency dropped -> inf volume) lived in all three copies;
// the fix had to be applied three times, and the stitch/mpi tie-break silently drifted (ENH-5).
//
// This centralizes the three pieces that must agree:
//   * the reduce/tie-break  -- keep the pair's LOWEST out_elev; break an elevation tie by the LOWER
//     out_cell (order-independent, so intra-then-seam visitation order can't change the result);
//   * the skip rule         -- skip a cell iff it is NoData AND not OCEAN (a NoData-as-OCEAN cell is
//     real base level, so its adjacency to land is a genuine basin->ocean outlet);
//   * the intra-tile D8 scan and the seam HandleEdge match.
// Call sites pass accessors (label/elev/cell-index/nodata, by grid position or by strip row), so the
// same code serves the global-grid callers (stitch, oracle) and the tile-local + Comm-strip caller
// (distributed). O(cells) + O(boundary).
#pragma once

#include <dephier/dephier.hpp>

#include <algorithm>
#include <map>
#include <utility>

namespace dh = richdem::dephier;

// The skip rule, defined once: a cell participates in an outlet only if it is not "NoData-and-not-OCEAN".
inline bool outlet_skip(bool is_nodata, dh::dh_label_t label){ return is_nodata && label!=dh::OCEAN; }

// One reduced outlet DB: for each depression pair, the outlet is the HIGHER cell of an adjacency; we keep
// the pair's LOWEST such out_elev, breaking an elevation tie by the LOWER out_cell. Templated on the cell-
// index type (the stitch uses dh::flat_c_idx = uint32_t; the mpi uses global int64_t indices).
template<class CellIdx>
struct OutletDB {
  using Key = std::pair<dh::dh_label_t, dh::dh_label_t>;
  std::map<Key, std::pair<float,CellIdx>> db;

  void reduce(dh::dh_label_t la, dh::dh_label_t lb, float ea, CellIdx ca, float eb, CellIdx cb){
    if(la==lb) return;                                     // same depression: not an outlet
    const float   oe = (ea>=eb)? ea : eb;                  // outlet = the higher of the pair
    const CellIdx oc = (ea>=eb)? ca : cb;
    const auto key = std::minmax(la, lb);
    const auto it = db.find({key.first,key.second});
    if(it==db.end() || oe < it->second.first || (oe==it->second.first && oc < it->second.second))
      db[{key.first,key.second}] = {oe, oc};
  }
};

// Intra-tile D8 scan over columns [xlo,xhi) x rows [0,H) on a WxH grid. `same_tile(x,nx)` gates each
// neighbour to the focal cell's tile (cross-seam pairs are handled by outlet_scan_seam). All accessors
// take GLOBAL grid coordinates (x,y); a tile-local caller maps them inside its lambdas.
template<class CellIdx, class Label,class Elev,class Cidx,class Nodata,class SameTile>
void outlet_scan_intra(OutletDB<CellIdx>& out, int xlo,int xhi,int W,int H,
                       Label label, Elev elev, Cidx cidx, Nodata nodata, SameTile same_tile){
  for(int y=0;y<H;y++) for(int x=xlo;x<xhi;x++){
    if(outlet_skip(nodata(x,y), label(x,y))) continue;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
      if(!dx && !dy) continue;
      const int nx=x+dx, ny=y+dy;
      if(nx<0||nx>=W||ny<0||ny>=H) continue;
      if(!same_tile(x,nx)) continue;                       // neighbour must be in the focal cell's tile
      if(outlet_skip(nodata(nx,ny), label(nx,ny))) continue;
      out.reduce(label(x,y),label(nx,ny), elev(x,y),cidx(x,y), elev(nx,ny),cidx(nx,ny));
    }
  }
}

// Seam HandleEdge match across ONE seam: side-A cell at row y vs side-B cells at rows y-1,y,y+1 (Barnes'
// parallel priority-flood edge pairing). Each side fixes a column and supplies label/elev/cell-index/
// nodata BY ROW -- serving both the grid-column pair (stitch, oracle) and own-column-vs-Comm-strip
// (distributed), whose side B is the neighbour's exchanged left-edge strip.
template<class CellIdx, class LA,class EA,class CA,class NA, class LB,class EB,class CB,class NB>
void outlet_scan_seam(OutletDB<CellIdx>& out, int H,
                      LA lblA,EA elvA,CA cidA,NA ndA,
                      LB lblB,EB elvB,CB cidB,NB ndB){
  for(int y=0;y<H;y++){
    if(outlet_skip(ndA(y), lblA(y))) continue;
    for(int ny=y-1;ny<=y+1;ny++){
      if(ny<0||ny>=H) continue;
      if(outlet_skip(ndB(ny), lblB(ny))) continue;
      out.reduce(lblA(y),lblB(ny), elvA(y),cidA(y), elvB(ny),cidB(ny));
    }
  }
}
