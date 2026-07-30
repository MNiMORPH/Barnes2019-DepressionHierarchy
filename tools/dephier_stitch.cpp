// dephier_stitch -- in-process tiled DepressionHierarchy build, diffed against serial.
//
// Stitch with BOUNDARY seeding (PARALLEL_DEPHIER_PLAN.md §3, §7.1):
//
//   per tile:  extract sub-DEM -> label ocean; pre-label seam-edge cells whose
//              steepest descent CROSSES the seam as BOUNDARY (so they seed as
//              provisional exterior and never become spurious pits) -> FloodAndAssignDepressions
//   global:    remap tiles into one namespace (shared ocean, offset the rest);
//              resolve every BOUNDARY cell to the real depression it drains into via a
//              boundary-graph pass (per-tile-local walk + cross-tile chaining, below);
//              build the outlet set as per-tile intra-tile outlets + Barnes' HandleEdge
//              across the seam perimeters; run one global ConstructHierarchyAndVolumes.
//   verify:    canonical signature must equal the serial build's.
//
// Both the outlet step and the conduit resolution are in the distributable shape: each
// tile works from its own arrays and only the 1-cell perimeter strips cross the seam
// (Barnes' HandleEdge; drain() at seam seeds). In this in-process harness that is
// organizationally equivalent to a single global pass, but no step reads another tile's
// interior -- so the per-rank footprint is O(N/P) + O(boundary). The one part still
// centralized is the final ConstructHierarchyAndVolumes (grid-free, O(#depressions)); a fully-distributed
// 2016-style join is future work (PARALLEL_DEPHIER_PLAN.md §10).

#include "dh_canonical.hpp"
#include "dh_collapse.hpp"   // CollapseSeamArtifacts       (shared with tools/dephier_mpi.cpp)
#include "dh_flats.hpp"      // ResolveFlatFlowdirs[_*]   (shared with tools/dephier_mpi.cpp)
#include "dh_outlets.hpp"    // OutletDB / outlet_scan_*    (shared with tools/dephier_mpi.cpp, ENH-5)

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>
#include <richdem/flowmet/d8_flowdirs.hpp>
#include <richdem/flats/flat_resolution.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;

using dh::dh_label_t;
using dh::OCEAN;
using dh::NO_DEP;
using dh::BOUNDARY;

static rd::Array2D<float> extract_cols(const rd::Array2D<float> &full, int x0, int x1){
  rd::Array2D<float> t(x1 - x0, full.height());
  t.setNoData(full.noData());
  for(int y=0;y<full.height();y++)
    for(int x=x0;x<x1;x++)
      t(x - x0, y) = full(x, y);
  return t;
}

static rd::Array2D<dh_label_t> ocean_labels(const rd::Array2D<float> &dem, float ocean_level){
  rd::Array2D<dh_label_t> label(dem.width(), dem.height(), NO_DEP);
  for(unsigned int i=0;i<dem.size();i++)
    if(dem.isNoData(i) || dem(i)==ocean_level)
      label(i) = OCEAN;
  return label;
}

template<class elev_t>
struct Tile {
  rd::Array2D<elev_t>             dem;
  rd::Array2D<dh_label_t>         label;
  rd::Array2D<int8_t>             fd;      // the tile's own flood flowdirs (tile-local)
  dh::DepressionHierarchy<elev_t> deps;
  std::vector<dh::Outlet<elev_t>> outlets;
  int        x0 = 0;
  dh_label_t offset = 0;
  dh_label_t g(dh_label_t local) const {
    if(local==OCEAN || local==BOUNDARY) return local;   // sentinels are global as-is
    return offset + (local - 1);
  }
};

// Whole-grid stitch state: the environment (shared, immutable) plus the built objects threaded
// through the build stages below. This mirrors dephier_mpi's RankCtx -- there each rank carries one
// tile's slice; here a single process holds the whole grid, so the "state" is the global tree plus
// the label/flowdir grids. The build stages read/write these members; main runs them in order and
// then reads the results (via reference-aliases) for the verification section.
struct StitchState {
  const rd::Array2D<float>& full;               // environment (immutable through the build)
  const std::vector<int>&   bounds;
  int   ntiles, W, H;
  float ocean_level;
  int   halo_cap;
  std::vector<Tile<float>>       tiles;         // built state, filled as the pipeline runs
  dh_label_t                     n_global = 1;  // 1 + total depressions across tiles
  dh::DepressionHierarchy<float> G;             // assembled global hierarchy
  rd::Array2D<dh_label_t>        gLabel;        // global label grid (resolved, then flat-reconciled)
  rd::Array2D<int8_t>            gFix;          // seam-fixed per-cell flowdirs
  std::vector<dh::Outlet<float>> outlets;       // re-derived outlet set fed to ConstructHierarchyAndVolumes

  StitchState(const rd::Array2D<float>& full_, const std::vector<int>& bounds_, float ocean_level_, int halo_cap_)
    : full(full_), bounds(bounds_), ntiles((int)bounds_.size()-1), W(full_.width()), H(full_.height()),
      ocean_level(ocean_level_), halo_cap(halo_cap_) {}

  int tile_of(int x) const { int t=0; while(t+1<(int)bounds.size() && x>=bounds[t+1]) t++; return t; }
  bool is_ocean(int x,int y) const { return full.isNoData(x,y) || full(x,y)==ocean_level; }
  // Lowest strictly-downhill LAND neighbour on the full grid, and whether the cell touches ocean.
  // An ocean neighbour is the sea (effectively -inf), so a cell touching ocean drains to the ocean
  // regardless of its land neighbours. Ties (equal-lowest neighbours) are broken by HIGHEST cell
  // index, matching the flood's radix pop order (radix_heap sorts each equal-elevation bucket
  // ascending and pops from the back, so the higher-index cell pops first and claims the shared
  // upslope neighbour). This matters at a seam: a divide cell with a tied descent -- one neighbour
  // in each tile -- is claimed by the serial flood from its highest-index neighbour, so drain() must
  // pick the same side or the cross-seam case is mislabelled.
  bool drain(int x,int y,int &lx,int &ly,bool &to_ocean) const {
    const float focal = full(x,y);
    float best = std::numeric_limits<float>::infinity(); bool found=false; to_ocean=false;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){        // index-ascending order
      if(!dx && !dy) continue;
      int nx=x+dx, ny=y+dy;
      if(!full.inGrid(nx,ny)) continue;
      if(is_ocean(nx,ny)){ to_ocean=true; continue; }
      const float e = full(nx,ny);
      if(e >= focal) continue;                                  // only strictly downhill
      if(e <= best){ best=e; lx=nx; ly=ny; found=true; }        // <= : later (higher-index) tie wins
    }
    return found;
  }
};

int main(int argc, char **argv){
  if(argc!=4 && argc!=5){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level> <Split Cols (comma-sep)> [flat halo cap]\n";
    return -1;
  }
  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);
  const int         halo_cap    = (argc==5) ? std::stoi(argv[4]) : std::numeric_limits<int>::max();

  rd::Array2D<float> full(in_name);

  std::vector<int> bounds = {0};
  {
    std::string s = argv[3];
    size_t p = 0;
    while(true){
      size_t q = s.find(',', p);
      bounds.push_back(std::stoi(s.substr(p, q-p)));
      if(q==std::string::npos) break;
      p = q+1;
    }
  }
  bounds.push_back(full.width());
  const int ntiles = bounds.size() - 1;

  const auto tile_of = [&](int x){ int t=0; while(t+1<(int)bounds.size() && x>=bounds[t+1]) t++; return t; };
  const auto is_ocean = [&](int x,int y){ return full.isNoData(x,y) || full(x,y)==ocean_level; };
  // Lowest strictly-downhill LAND neighbour on the full grid, and whether the cell
  // touches ocean. An ocean neighbour is the sea (effectively -inf), so a cell
  // touching ocean drains to the ocean regardless of its land neighbours.
  //
  // Ties (equal-lowest neighbours) are broken by HIGHEST cell index, matching the
  // flood's radix pop order (radix_heap sorts each equal-elevation bucket ascending
  // and pops from the back, so the higher-index cell pops first and claims the shared
  // upslope neighbour). This matters at a seam: a divide cell with a tied descent —
  // one neighbour in each tile — is claimed by the serial flood from its highest-index
  // neighbour, so drain() must pick the same side or the cross-seam case is mislabelled.
  const auto drain = [&](int x,int y,int &lx,int &ly,bool &to_ocean)->bool{
    const float focal = full(x,y);
    float best = std::numeric_limits<float>::infinity(); bool found=false; to_ocean=false;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){        // index-ascending order
      if(!dx && !dy) continue;
      int nx=x+dx, ny=y+dy;
      if(!full.inGrid(nx,ny)) continue;
      if(is_ocean(nx,ny)){ to_ocean=true; continue; }
      const float e = full(nx,ny);
      if(e >= focal) continue;                                  // only strictly downhill
      if(e <= best){ best=e; lx=nx; ly=ny; found=true; }        // <= : later (higher-index) tie wins
    }
    return found;
  };

  // ---- per-tile: pre-label seam-crossing edges BOUNDARY, then FloodAndAssignDepressions ----
  // Each tile keeps its OWN flowdirs (tile.fd) -- the flood's drainage map, tile-local.
  // Conduit resolution below follows each tile's own fd (no global flow map), so it is
  // internally consistent with the flood -- including flats, where the flood's
  // deterministic order (radix_heap bucket sort) fixes the direction.
  std::vector<Tile<float>> tiles(ntiles);
  dh_label_t next_offset = 1;
  for(int t=0;t<ntiles;t++){
    auto &tile = tiles[t];
    tile.x0    = bounds[t];
    tile.dem   = extract_cols(full, bounds[t], bounds[t+1]);
    tile.label = ocean_labels(tile.dem, ocean_level);

    // internal-seam edge columns of this tile
    std::vector<int> seam_cols;
    if(t>0)          seam_cols.push_back(bounds[t]);       // left edge
    if(t<ntiles-1)   seam_cols.push_back(bounds[t+1]-1);   // right edge
    for(int gc : seam_cols)
      for(int y=0;y<full.height();y++){
        const int lc = gc - tile.x0;
        if(tile.label(lc,y)!=NO_DEP) continue;             // ocean stays ocean
        int lx,ly; bool to_ocean;
        // BOUNDARY only if it drains across the seam and NOT to the ocean
        if(drain(gc,y,lx,ly,to_ocean) && !to_ocean && tile_of(lx)!=t)
          tile.label(lc,y) = BOUNDARY;
      }

    tile.fd = rd::Array2D<int8_t>(tile.dem.width(), tile.dem.height(), rd::NO_FLOW);
    // A tile may be a bowl interior (no base-level seed); permit the pit-only flood. Its open
    // top depression is closed later across the seam by HandleEdge + ConstructHierarchyAndVolumes (ENH-2).
    dh::FloodAndAssignDepressions<float,rd::Topology::D8>(tile.dem, tile.label, tile.fd, tile.deps, tile.outlets, /*permit_without_baselevel_seed=*/true);
    tile.offset = next_offset;
    next_offset += tile.deps.size() - 1;
  }
  const dh_label_t n_global = next_offset;

  // ---- assemble global depressions + global label grid ----
  dh::DepressionHierarchy<float> G(n_global);
  G[0] = tiles[0].deps[0];
  G[0].dep_label = 0;
  for(auto &tile : tiles)
    for(dh_label_t k=1;k<tile.deps.size();k++){
      auto d = tile.deps[k];
      d.dep_label = tile.g(k);
      if(d.pit_cell!=dh::NO_VALUE){
        int lx,ly; tile.dem.iToxy(d.pit_cell, lx, ly);
        d.pit_cell = full.xyToI(tile.x0 + lx, ly);
      }
      G[tile.g(k)] = d;
    }
  rd::Array2D<dh_label_t> gLabel(full.width(), full.height(), OCEAN);
  for(auto &tile : tiles)
    for(int y=0;y<tile.dem.height();y++)
      for(int x=0;x<tile.dem.width();x++)
        gLabel(tile.x0 + x, y) = tile.g(tile.label(x, y));

  // ---- conduit resolution: distributable boundary-graph pass (PARALLEL_DEPHIER_PLAN.md §3) ----
  // A BOUNDARY cell is a through-flowing cell whose drainage leaves its tile; we must
  // find the real depression (or ocean) it ends in WITHOUT walking the full grid.
  //   Phase 1 (per tile, LOCAL): follow the tile's OWN flowdirs from each BOUNDARY cell
  //     to a tile-local terminal (a depression/ocean -> its global label) or a seam exit
  //     (the cross-seam cell it hands off to). Touches only tiles[t].label / tiles[t].fd
  //     and the 1-cell cross-seam neighbour via drain() -- the perimeter strip a rank
  //     already exchanges. A tile's own fd never points across a seam (FloodAndAssignDepressions sees only
  //     in-tile neighbours), so the walk stays in-tile until it terminates or hits a
  //     NO_FLOW seam seed, where it crosses.
  //   Phase 2 (CHAIN): follow those exits across tiles until a terminal. Conduit paths
  //     are shallow (measured <= 6 tile-crossings even at 8 tiles), so this converges in
  //     a few hops, reading only the O(boundary) Phase-1 results plus edge labels --
  //     never another tile's interior. Order-independent: Phase-1 walks the unmodified
  //     tile labels; the resolved global labels are written into gLabel together at the
  //     end. This is the footprint-bounded form of what a single full-grid walk did.
  {
    struct LW { bool terminal; dh_label_t label; int ex, ey; };   // terminal label, else exit target (global x,y)
    const auto localwalk = [&](int t, int gx, int gy)->LW {
      const auto &tile = tiles[t];
      for(unsigned guard=0; guard<=tile.label.size(); guard++){
        const dh_label_t ll = tile.label(gx - tile.x0, gy);
        if(ll!=BOUNDARY) return { true, tile.g(ll), 0, 0 };         // local depression or ocean
        const int fd = tile.fd(gx - tile.x0, gy);
        if(fd!=rd::NO_FLOW){ gx += rd::d8x[fd]; gy += rd::d8y[fd]; } // stays within the tile
        else {                                                      // seam seed: cross the strip
          int nx,ny; bool to_ocean;
          if(!drain(gx,gy,nx,ny,to_ocean) || to_ocean) return { true, OCEAN, 0, 0 };
          return { false, OCEAN, nx, ny };                          // hand off to the neighbour tile
        }
      }
      return { true, OCEAN, 0, 0 };                                 // safety net (no cycle expected)
    };

    // Phase 1: each tile resolves its own BOUNDARY cells with a tile-local walk.
    std::map<std::pair<int,int>, LW> res;
    for(int t=0;t<ntiles;t++)
      for(int y=0;y<full.height();y++)
        for(int gx=tiles[t].x0; gx<bounds[t+1]; gx++)
          if(tiles[t].label(gx - tiles[t].x0, y)==BOUNDARY)
            res[{gx,y}] = localwalk(t, gx, y);

    // Phase 2: chain exits through the published results until a terminal.
    const auto chase = [&](int gx, int gy)->dh_label_t {
      for(unsigned guard=0; guard<=res.size()+1; guard++){
        const LW &r = res.at({gx,gy});
        if(r.terminal) return r.label;
        const int nt = tile_of(r.ex);
        const dh_label_t nl = tiles[nt].label(r.ex - tiles[nt].x0, r.ey);
        if(nl!=BOUNDARY) return tiles[nt].g(nl);   // entry cell is itself a terminal
        gx = r.ex; gy = r.ey;                       // else it is a BOUNDARY cell -> chain on its result
      }
      return OCEAN;
    };
    std::vector<std::pair<int,int>> bcells;
    bcells.reserve(res.size());
    for(const auto &kv : res) bcells.push_back(kv.first);
    std::vector<dh_label_t> resolved(bcells.size());
    for(size_t i=0;i<bcells.size();i++) resolved[i] = chase(bcells[i].first, bcells[i].second);
    for(size_t i=0;i<bcells.size();i++) gLabel(bcells[i].first, bcells[i].second) = resolved[i];
  }

  // ---- flowdir fix-up: restore serial-identical per-cell flowdirs at seam crossings ----
  // A tile flood cannot point a cell across its own boundary, so a cell whose true
  // drainage crosses the seam is left with a tile-local flowdir (NO_FLOW at a BOUNDARY
  // seam seed). Every other cell already matches serial (measured: divergence is 100%
  // seam-confined). The conduit pass already found each seam seed's cross-seam drain
  // target, so point the seed's flowdir there -- the same steepest-descent, highest-index-
  // tie direction the serial flood assigns. This is an O(boundary) edit using only the
  // perimeter strip; here it produces a global field for validation against serial.
  rd::Array2D<int8_t> gFix(full.width(), full.height(), rd::NO_FLOW);
  for(auto &tile : tiles)
    for(int y=0;y<tile.dem.height();y++)
      for(int x=0;x<tile.dem.width();x++)
        gFix(tile.x0 + x, y) = tile.fd(x, y);
  const auto dir_to = [&](int dx,int dy)->int8_t{                 // D8 step -> flowdir index
    for(int n=1;n<=8;n++) if(rd::d8x[n]==dx && rd::d8y[n]==dy) return (int8_t)n;
    return rd::NO_FLOW;
  };
  for(auto &tile : tiles)
    for(int y=0;y<tile.dem.height();y++)
      for(int x=0;x<tile.dem.width();x++){
        const int gx=tile.x0+x;
        if(tile.fd(x,y)==rd::NO_FLOW && !full.isNoData(gx,y) && tile.label(x,y)!=OCEAN){
          // A cell the tile flood left as a local pit (NO_FLOW) but whose true drainage
          // exits across the seam: a land seam seed (strictly-lower cross neighbour) or a
          // seam-straddling flat (equal cross neighbour). Serial claims it from its lowest
          // cross-tile neighbour AT OR BELOW its elevation, ties by highest cell index
          // (matching the flood). A genuine global pit has no such neighbour -> stays
          // NO_FLOW, as in serial.
          const float fe = full(gx,y); const long fi = full.xyToI(gx,y);
          int bnx=0,bny=0; long bi=-1; float be=std::numeric_limits<float>::infinity();
          for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            if(!dx && !dy) continue;
            const int nx=gx+dx, ny=y+dy;
            if(!full.inGrid(nx,ny) || full.isNoData(nx,ny)) continue;
            if(tile_of(nx)==tile_of(gx)) continue;                     // cross-seam neighbours only
            const float e=full(nx,ny); if(e>fe) continue;             // at or below the pit
            const long idx=full.xyToI(nx,ny);
            if(e<be || (e==be && idx>bi)){ be=e; bi=idx; bnx=nx; bny=ny; }
          }
          // Apply if that neighbour is strictly lower, or an equal-elevation cell of
          // HIGHER index -- the flood makes the highest-index cell of a straddling flat
          // the pit and points the others at it, so an equal lower-index neighbour means
          // THIS cell is the pit (leave it NO_FLOW).
          if(bi>=0 && (be<fe || bi>fi)) gFix(gx,y) = dir_to(bnx-gx, bny-y);
        } else if(tile.label(x,y)==OCEAN && !full.isNoData(gx,y)){    // sea-draining cell
          // serial claims a sea-draining cell from its HIGHEST-index adjacent ocean cell
          // (all ocean seeds, popped highest-index-first). When that lies across the seam
          // the tile picks a different one; restore serial's choice.
          int bi=-1, bnx=0, bny=0;
          for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            if(!dx && !dy) continue;
            const int nx=gx+dx, ny=y+dy;
            if(!full.inGrid(nx,ny) || !is_ocean(nx,ny)) continue;
            const int idx=full.xyToI(nx,ny);
            if(idx>bi){ bi=idx; bnx=nx; bny=ny; }
          }
          if(bi>=0 && tile_of(bnx)!=tile_of(gx))                       // highest ocean nbr is across seam
            gFix(gx,y) = dir_to(bnx-gx, bny-y);
        }
      }
  // Flat cells: replace the flood's order-dependent claim with Barnes-2014 resolve_flats
  // (deterministic from geometry -> agrees with serial). MOVE 2: footprint-bounded, per-tile
  // resolve_flats with an adaptive boundary halo (no full grid); bit-identical to full-grid.
  const int flat_capped = ResolveFlatFlowdirsAdaptiveHalo(full, bounds, gFix, halo_cap);
  if(flat_capped) std::cerr<<"flat resolution: "<<flat_capped<<" tile(s) hit the halo cap ("<<halo_cap
                           <<") -- valid but possibly not serial-identical in giant-flat interiors\n";

  // ---- global outlet set, in the distributable shape (Barnes' join): each tile
  // derives its own outlets from its resolved labels (intra-tile adjacencies), and
  // only the perimeter strips cross the seam via HandleEdge. Both feed one database
  // keyed on the depression pair, keeping the lowest max-of-pair -- the serial PhaseB
  // rule. Intra-tile + cross-seam together cover every adjacency, so the result is
  // identical to a single global pass. ----
  //
  // WHY re-derive here instead of reusing FloodAndAssignDepressions's per-tile `tile.outlets` (this is NOT fluffy
  // duplication -- it is load-bearing for the distributed build):
  //   * FloodAndAssignDepressions found its outlets in each tile's PRE-conduit labels (BOUNDARY sentinels included).
  //     The conduit pass then resolves every BOUNDARY cell to its true drain target, which changes
  //     *which* depressions actually meet. Re-deriving from the RESOLVED labels captures that for
  //     free; remapping `tile.outlets` through the resolution would be fiddly and error-prone.
  //   * In the real MPI build each rank re-derives from ITS OWN resolved labels locally (dephier_mpi
  //     does the same from glab_pc) -- no raw outlet set is shipped or merged. That locality is the
  //     whole point; `tile.outlets` is deliberately dropped.
  // THE COST (record it, so the risk is visible): this is a SECOND outlet-discovery path, and it must
  // reproduce FloodAndAssignDepressions's semantics faithfully or the two drift. One drift was a real bug: a NoData cell
  // is OCEAN (ocean_labels), so a basin->NoData-ocean adjacency is a genuine outlet; the scan used to
  // `continue` on all NoData cells and drop it, leaving the basin unclosed (inf volume). The `skip()`
  // below closes that gap (skip a cell only if NoData AND not OCEAN). A debug flag audits for further
  // drift by diffing this set against FloodAndAssignDepressions's tile.outlets -- see DH_AUDIT_OUTLETS below.
  // ---- FLAT-PARTITION REPLAY: reconstruct serial's exact cell->leaf partition (flag-gated; ROAD A) ----
  // SPLIT_INVARIANT_FLATS_PLAN.md ROAD A. The tiled flood labels a divide/flat cell by a SEAM-DEPENDENT
  // wavefront, so a cell can land in the wrong leaf -> wrong cell_count/dep_vol (the cell-assignment DIFFER
  // class). Serial's label = the depression the FLOOD claims the cell into, which for flats is a geodesic
  // partition tie-broken by the radix pop order -- NOT the pit a flowdir drains to (that mistake emptied
  // sill-flat leaves; see the containment cc-pass finding). Here we REPLAY serial's flood-labelling exactly
  // (proven bit-identical partition on 49 DEMs incl. Corsica + adversarial fractals; tools/
  // flat_partition_replay_proof.cpp): seeds = every OCEAN cell + every land cell with NO strictly-lower
  // neighbour (the land_seed/pit set, dephier.hpp:412 -- includes all flat interiors), each at its dem
  // elevation; process elevation buckets ASCending, within a bucket HIGHEST-index first with same-elevation
  // labels appended and popped LIFO (the fork radix_heap tie-break, radix_heap.hpp:391); a popped NO_DEP cell
  // starts a new pit; ocean at its dem elevation, propagating OCEAN. Then map each replay-basin onto G's
  // EXISTING leaf via that leaf's pit cell (the min leaf when a seam split one basin into several -> they
  // merge) and overwrite gLabel; compact G. Structure follows: outlets are re-derived
  // from this corrected gLabel, so ConstructHierarchyAndVolumes builds on serial's partition. Full-grid here (the flag's
  // reproducibility trade); distributable later via ENH-1 seam-exchange replaying the flood ORDER across the
  // seam. Default OFF => byte-identical.
  const bool flat_replay = std::getenv("DH_FLAT_PARTITION_REPLAY")!=nullptr;
  if(flat_replay){
    const int W=full.width(), H=full.height();
    // 1. Full pit-up replay of the flood's label partition into r (fresh namespace; OCEAN kept as OCEAN).
    rd::Array2D<dh_label_t> r(W,H,dh::NO_DEP);
    std::map<float,std::vector<long>> buckets;
    const auto rpush=[&](int x,int y){ buckets[full(x,y)].push_back(full.xyToI(x,y)); };
    for(int y=0;y<H;y++) for(int x=0;x<W;x++){
      if(is_ocean(x,y)){ r(x,y)=OCEAN; rpush(x,y); continue; }
      const float e=full(x,y); bool has_lower=false;
      for(int dy=-1;dy<=1&&!has_lower;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=x+dx,ny=y+dy;
        if(full.inGrid(nx,ny)&&full(nx,ny)<e){ has_lower=true; break; } }
      if(!has_lower) rpush(x,y);
    }
    dh_label_t nextr=1;
    std::vector<long> rb_pit(1, -1);                          // rb_pit[label] = the cell that BORE that basin
    while(!buckets.empty()){
      auto it=buckets.begin(); const float e=it->first;
      std::vector<long> cur=std::move(it->second); buckets.erase(it);
      std::sort(cur.begin(),cur.end());                       // ascending index; pop_back = highest first
      while(!cur.empty()){
        const long ci=cur.back(); cur.pop_back();
        int cx,cy; full.iToxy(ci,cx,cy);
        dh_label_t cl=r(cx,cy);
        if(cl==dh::NO_DEP){ cl=nextr++; r(cx,cy)=cl; rb_pit.push_back(ci); }   // popped NO_DEP -> new pit (= serial's pit)
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=cx+dx,ny=cy+dy;
          if(!full.inGrid(nx,ny) || r(nx,ny)!=dh::NO_DEP) continue;
          r(nx,ny)=cl; if(full(nx,ny)==e) cur.push_back(full.xyToI(nx,ny)); else buckets[full(nx,ny)].push_back(full.xyToI(nx,ny)); }
      }
    }
    // 2. Map each replay-basin to G's EXISTING leaf via that leaf's pit cell (min leaf = canonical; several
    //    leaves in one replay-basin = a seam-split basin -> they merge into the min).
    std::vector<dh_label_t> rb2leaf(nextr, dh::NO_VALUE);
    for(dh_label_t i=1;i<G.size();i++){
      if(G[i].pit_cell==dh::NO_VALUE) continue;
      int px,py; full.iToxy(G[i].pit_cell,px,py); const dh_label_t rb=r(px,py);
      if(rb==OCEAN || rb>=nextr) continue;
      if(rb2leaf[rb]==dh::NO_VALUE || i<rb2leaf[rb]) rb2leaf[rb]=i;
    }
    // Stamp each canonical leaf with the replay's own pit (= serial's pit: the lowest, highest-index cell of
    // the basin). pit_elev IS in the canonical signature, so a merged leaf must carry the basin's true min,
    // not the arbitrary min-label fragment's pit (fixes testdem8 sp3: serial pit_elev 2 vs fragment's 4).
    for(dh_label_t rb=1;rb<nextr;rb++){
      const dh_label_t L=rb2leaf[rb];
      if(L==dh::NO_VALUE) continue;
      G[L].pit_cell=(dh::flat_c_idx)rb_pit[rb]; int px,py; full.iToxy(rb_pit[rb],px,py); G[L].pit_elev=full(px,py);
    }
    // 3. Overwrite gLabel with serial's partition (in G's namespace).
    long relabeled=0, orphan=0;
    for(int y=0;y<H;y++) for(int x=0;x<W;x++){
      if(is_ocean(x,y)){ gLabel(x,y)=OCEAN; continue; }
      const dh_label_t rb=r(x,y);
      if(rb==OCEAN){ if(gLabel(x,y)!=OCEAN){ gLabel(x,y)=OCEAN; relabeled++; } continue; }
      const dh_label_t L = (rb<nextr) ? rb2leaf[rb] : dh::NO_VALUE;
      if(L==dh::NO_VALUE){ orphan++; continue; }               // no G-leaf for this basin: keep existing label
      if(gLabel(x,y)!=L){ gLabel(x,y)=L; relabeled++; }
    }
    // 4. Compact G: keep only leaves still referenced by gLabel (dropped = merged-away / emptied), renumber.
    std::vector<char> used(G.size(),0); used[0]=1;             // ocean node always kept
    for(int y=0;y<H;y++) for(int x=0;x<W;x++) if(gLabel(x,y)!=OCEAN) used[gLabel(x,y)]=1;
    std::vector<dh_label_t> dense(G.size(),0); dh_label_t nn=0;
    for(dh_label_t i=0;i<G.size();i++) if(used[i]) dense[i]=nn++;
    dh::DepressionHierarchy<float> G2(nn);
    for(dh_label_t i=0;i<G.size();i++) if(used[i]){ G2[dense[i]]=G[i]; G2[dense[i]].dep_label=dense[i]; }
    for(int y=0;y<H;y++) for(int x=0;x<W;x++) if(gLabel(x,y)!=OCEAN) gLabel(x,y)=dense[gLabel(x,y)];
    const dh_label_t old_leaves=G.size(); G = std::move(G2);
    std::cerr<<"FLAT-PARTITION-REPLAY: "<<nextr-1<<" replay basins, relabeled "<<relabeled
             <<" cells, leaves "<<old_leaves<<"->"<<nn<<" orphan_cells="<<orphan<<"\n";
  }

  std::vector<dh::Outlet<float>> outlets;
  {
    // Re-derive the outlet set from the resolved label grid via the shared scan (dh_outlets.hpp, ENH-5):
    // intra-tile D8 adjacencies, then Barnes' HandleEdge across each internal seam. The reduce keeps the
    // pair's lowest out_elev (tie -> lower out_cell), and OutletSkip drops a cell only if NoData-and-not-
    // OCEAN -- the two rules that must match the mpi oracle/distributed scans (and once did not).
    const int W = full.width(), H = full.height();
    OutletDB<dh::flat_c_idx> odb;
    const auto label  = [&](int x,int y){ return gLabel(x,y); };
    const auto elev   = [&](int x,int y){ return full(x,y); };
    const auto cidx   = [&](int x,int y){ return full.xyToI(x,y); };
    const auto nodata = [&](int x,int y){ return full.isNoData(x,y); };
    OutletScanIntra(odb, 0, W, W, H, label, elev, cidx, nodata,
                      [&](int x,int nx){ return tile_of(x)==tile_of(nx); });
    for(size_t b=1;b+1<bounds.size();b++){
      const int cA=bounds[b]-1, cB=bounds[b];     // the two columns straddling the seam
      OutletScanSeam(odb, H,
        [&](int y){return gLabel(cA,y);},[&](int y){return full(cA,y);},[&](int y){return full.xyToI(cA,y);},[&](int y){return full.isNoData(cA,y);},
        [&](int y){return gLabel(cB,y);},[&](int y){return full(cB,y);},[&](int y){return full.xyToI(cB,y);},[&](int y){return full.isNoData(cB,y);});
    }

    for(auto &kv : odb.db){
      dh::Outlet<float> o;
      o.depa = kv.first.first; o.depb = kv.first.second;
      o.out_elev = kv.second.first; o.out_cell = kv.second.second;
      outlets.push_back(o);
    }
  }

  // ---- DEBUG: audit the re-derived global outlet set against FloodAndAssignDepressions's own per-tile outlets ----
  // Set DH_AUDIT_OUTLETS=1 to diff the two outlet-discovery paths in real time. The re-derivation
  // (above) must reproduce every outlet FloodAndAssignDepressions found, or the two drift (that drift is how the
  // NoData-ocean inf bug arose). For each tile we remap FloodAndAssignDepressions's tile.outlets to global labels via
  // tile.g and check the pair appears in the re-derived set with the same outlet elevation. Pairs
  // touching BOUNDARY are skipped -- their endpoints are only meaningful AFTER conduit resolution,
  // which is exactly what the re-derivation applies, so a raw comparison there is not apples-to-apples.
  if(std::getenv("DH_AUDIT_OUTLETS")){
    std::map<std::pair<dh_label_t,dh_label_t>,float> grid;
    for(const auto &o : outlets) grid[std::minmax(o.depa,o.depb)] = o.out_elev;
    long miss=0, ediff=0;
    for(const auto &tile : tiles)
      for(const auto &o : tile.outlets){
        const dh_label_t a=tile.g(o.depa), b=tile.g(o.depb);
        if(a==BOUNDARY || b==BOUNDARY || a==b) continue;         // BOUNDARY: needs resolution to compare
        const auto it = grid.find(std::minmax(a,b));
        if(it==grid.end()){ miss++;
          std::cerr<<"AUDIT miss: FloodAndAssignDepressions outlet {"<<a<<","<<b<<"} elev="<<o.out_elev<<" absent from re-derived set\n"; }
        else if(it->second!=o.out_elev){ ediff++;
          std::cerr<<"AUDIT elev: pair {"<<a<<","<<b<<"} FloodAndAssignDepressions="<<o.out_elev<<" re-derived="<<it->second<<"\n"; }
      }
    std::cerr<<"AUDIT-OUTLETS "<<in_name<<" splits="<<argv[3]<<": "<<miss<<" missing pair(s), "
             <<ediff<<" elevation diff(s) (FloodAndAssignDepressions tile.outlets vs re-derived global)\n";
  }

  // ---- DEBUG: compare the tiled outlet set fed to ConstructHierarchyAndVolumes against SERIAL's own outlets ----
  // Set DH_AUDIT_VS_SERIAL=1. ConstructHierarchy's outlet sort (dephier.hpp:660) is fully deterministic on the outlet
  // SET (out_elev -> out_cell -> endpoint pit cells), so if the tiled build feeds the SAME outlets serial
  // does, it reproduces serial's tree exactly. This audits that: a fresh serial FloodAndAssignDepressions's outlets vs the
  // tiled `outlets`, keyed namespace-free by the endpoint depressions' pit cells. A differing out_cell on a
  // shared pair (same out_elev) is a TIE the tiled outlet re-derivation broke differently than serial's
  // flood -- the source of the meta-ordering / meta-vs-ocean_linked residual, and an OUR-SIDE fix (no serial
  // change) iff serial's own choice is itself deterministic/reproducible.
  if(std::getenv("DH_AUDIT_VS_SERIAL")){
    rd::Array2D<dh_label_t> a_label = ocean_labels(full, ocean_level);
    rd::Array2D<int8_t> a_fd(full.width(), full.height(), rd::NO_FLOW);
    dh::DepressionHierarchy<float> a_deps; std::vector<dh::Outlet<float>> a_out;
    dh::FloodAndAssignDepressions<float,rd::Topology::D8>(full, a_label, a_fd, a_deps, a_out, true);
    const auto keyOf=[&](const dh::DepressionHierarchy<float>&D, dh_label_t a, dh_label_t b){
      const auto pc=[&](dh_label_t x)->long{ return (x==OCEAN||x>=D.size()||D[x].pit_cell==dh::NO_VALUE)?-1:(long)D[x].pit_cell; };
      const long pa=pc(a), pb=pc(b); return std::make_pair(std::min(pa,pb), std::max(pa,pb)); };
    std::map<std::pair<long,long>,std::pair<float,long>> smap, tmap;
    for(auto&o:a_out)   smap[keyOf(a_deps,o.depa,o.depb)]={o.out_elev,(long)o.out_cell};
    for(auto&o:outlets) tmap[keyOf(G,     o.depa,o.depb)]={o.out_elev,(long)o.out_cell};
    long same=0, cdiff=0, ediff=0, only_s=0, only_t=0;
    for(auto&kv:smap){ auto it=tmap.find(kv.first);
      if(it==tmap.end()){ only_s++; continue; }
      if(kv.second.first!=it->second.first) ediff++;
      else if(kv.second.second!=it->second.second){ cdiff++;
        std::cerr<<"  OUTLET out_cell diff at pit-pair ("<<kv.first.first<<","<<kv.first.second<<") elev="
                 <<kv.second.first<<": serial out_cell="<<kv.second.second<<" tiled="<<it->second.second<<"\n"; }
      else same++; }
    for(auto&kv:tmap) if(smap.find(kv.first)==smap.end()) only_t++;
    std::cerr<<"AUDIT-VS-SERIAL "<<in_name<<" splits="<<argv[3]<<": serial_outlets="<<a_out.size()
             <<" tiled_outlets="<<outlets.size()<<" same="<<same<<" out_cell_diff="<<cdiff
             <<" out_elev_diff="<<ediff<<" only_serial="<<only_s<<" only_tiled="<<only_t<<"\n";
  }

  // ---- one global ConstructHierarchyAndVolumes ----
  dh::ConstructHierarchyAndVolumes<float>(G, outlets, full, gLabel);

  // ---- §3.2 collapse pass: contract seam-split artifacts to a serial-identical tree ----
  // RETIRED when the flat-partition replay is on: the replay gives serial's exact partition, so outlets ->
  // ConstructHierarchyAndVolumes build serial's tree with no seam artifacts, and the compact already drops spurious leaves. The
  // collapse's meta-over-halves passes (seam-dependent) would then DISSOLVE REAL structure (measured:
  // kerry_test9 sp7 merged two genuine basins). This is the retirement the warning has been flagging.
  int n_collapsed = 0;
  if(!flat_replay){
    n_collapsed = CollapseSeamArtifacts(G, full, bounds);
    if(n_collapsed) std::cerr<<"collapse: contracted "<<n_collapsed<<" seam artifact(s)\n";
  }

  // ---- serial ground truth ----
  auto s_label = ocean_labels(full, ocean_level);
  rd::Array2D<int8_t> s_fd(full.width(), full.height(), rd::NO_FLOW);
  auto S = dh::GetDepressionHierarchy<float,rd::Topology::D8>(full, s_label, s_fd);
  // Serial reference uses the same deterministic flat routing (the proposed DH change:
  // resolve_flats instead of the flood's flat byproduct). The tree is untouched.
  ResolveFlatFlowdirs(full, s_fd);

  // ---- flowdir check: the fixed-up tiled flowdirs must equal serial's, cell for cell ----
  bool flowdir_ok = true;
  {
    long land=0, fdiff=0;
    for(unsigned int i=0;i<full.size();i++){
      if(full.isNoData(i)) continue;
      land++;
      if(gFix(i)!=s_fd(i)) fdiff++;
    }
    flowdir_ok = (fdiff==0);
    std::cout<<(flowdir_ok ? "FLOWDIR-MATCH " : "FLOWDIR-DIFFER ")<<in_name
             <<" fd_diff="<<fdiff<<"/"<<land<<"\n";
  }

  // ---- compare ----
  const std::string sig_stitch = dhtest::canonicalize(G);
  const std::string sig_serial = dhtest::canonicalize(S);
  const auto iv_stitch = dhtest::invariants(G);
  const auto iv_serial = dhtest::invariants(S);

  std::cerr<<"serial nodes="<<iv_serial.n_nodes<<" (leaf "<<iv_serial.n_leaf<<", meta "<<iv_serial.n_meta<<")"
           <<"  stitch nodes="<<iv_stitch.n_nodes<<" (leaf "<<iv_stitch.n_leaf<<", meta "<<iv_stitch.n_meta<<")\n";
  std::cerr<<"serial total_dep_vol="<<iv_serial.total_dep_vol<<"  stitch total_dep_vol="<<iv_stitch.total_dep_vol<<"\n";

  const bool ok = (sig_stitch==sig_serial);
  std::cout<<(ok ? "STITCH-MATCH " : "STITCH-DIFFER ")<<in_name<<" splits="<<argv[3]<<"\n";

  // Decomposition-correctness diagnostic: the NODE COUNT. The depression tree is split-invariant -- a
  // correct build has the same depressions no matter where the seams fall -- so an unequal node count is a
  // NECESSARY-condition symptom that the split changed the tree. It flags two things the diff must separate:
  // (1) a genuine seam artifact the collapse missed (a spurious extra depression -- a real bug to chase);
  // (2) the STRUCTURAL ConstructHierarchyAndVolumes tie-break sub-class -- at a TIED outlet, two basins rebuilt as a META vs one
  // ocean_linked into the other (same depressions + volume, but node count differs; kerry_test2 4 vs 3).
  // (2) changes serial output -> Richard-coordinated, acceptable under the volume-correct+valid-tree bar.
  // Note node count does NOT isolate (1) from all tie-break noise: the pure-ORDERING tie-break sub-class
  // preserves node count (DIFFER but DECOMP-CORRECT, e.g. kerry_test4 splits 3/8); only meta-vs-ocean_linked
  // moves it.
  const bool decomp_ok = (iv_stitch.n_nodes == iv_serial.n_nodes);
  std::cout<<(decomp_ok ? "STITCH-DECOMP-CORRECT " : "STITCH-DECOMP-INCORRECT ")<<in_name
           <<" splits="<<argv[3]<<" nodes(serial="<<iv_serial.n_nodes<<" stitch="<<iv_stitch.n_nodes<<")\n";

  // Volume-correctness verdict (VOL-MATCH/VOL-DIFFER): total depression volume vs serial. This is the
  // acceptance bar for cases that are NOT bit-identical to serial -- the bowl-interior (ENH-2) splits and
  // the tie-break class are volume-correct + valid-tree by design even when the tree signature DIFFERs, so
  // STITCH-MATCH can't assert them but this can. A small relative tolerance absorbs FP summation-order
  // differences (the distributed PhaseD sums per-rank partials in a different order than serial's single
  // pass); an open depression (total_dep_vol=inf) or a lost/merged basin makes the gap non-trivial -> DIFFER.
  const double v_s = iv_serial.total_dep_vol, v_d = iv_stitch.total_dep_vol;
  const bool vol_ok = std::isfinite(v_s) && std::isfinite(v_d) &&
                      std::fabs(v_d - v_s) <= 1e-6 * std::max(1.0, std::fabs(v_s));
  std::cout<<(vol_ok ? "STITCH-VOL-MATCH " : "STITCH-VOL-DIFFER ")<<in_name
           <<" splits="<<argv[3]<<" total_dep_vol(serial="<<v_s<<" stitch="<<v_d<<")\n";

  // Leaf-set verdict (LEAFSET-MATCH/DIFFER): does the tiled build have EXACTLY serial's multiset of LEAF
  // depressions -- same (pit_elev, out_elev, cell_count, dep_vol) -- ignoring the meta tree above them? This
  // separates the residual DIFFER classes from a real bug: if the leaf sets MATCH but STITCH-DIFFERs, the
  // ONLY difference is meta-tree SHAPE (the outlet-ordering + meta-vs-ocean_linked tie-break, both genuine
  // ConstructHierarchy ties). A LEAFSET-DIFFER on a STITCH-DIFFER case would be an actual cell-assignment/basin residual.
  {
    const auto leafset=[&](const dh::DepressionHierarchy<float>&D){
      std::vector<std::string> v;
      for(dh_label_t i=1;i<D.size();i++){ const auto&d=D[i];
        if(d.lchild==dh::NO_VALUE && d.rchild==dh::NO_VALUE)
          v.push_back(dhtest::quant(d.pit_elev,4)+","+dhtest::quant(d.out_elev,4)+","
                      +std::to_string(d.cell_count)+","+dhtest::quant(d.dep_vol,4)); }
      std::sort(v.begin(),v.end()); return v; };
    const auto ls=leafset(S), lg=leafset(G);
    std::cout<<((ls==lg)?"STITCH-LEAFSET-MATCH ":"STITCH-LEAFSET-DIFFER ")<<in_name
             <<" splits="<<argv[3]<<" leaves(serial="<<ls.size()<<" stitch="<<lg.size()<<")\n";
  }

  // Split-invariance signature: a hash of the tiled build's own canonical tree, printed ALWAYS. Two runs
  // at DIFFERENT splits must print the same STITCH-SIG iff the build is split-invariant on this DEM (the
  // real acceptance metric -- "identical DH regardless of how we split" -- which does NOT reference serial;
  // a valid tie-break vs serial keeps STITCH-MATCH DIFFER but leaves STITCH-SIG invariant across tilings).
  std::cout<<"STITCH-SIG "<<in_name<<" splits="<<argv[3]<<" sig="<<std::hash<std::string>{}(sig_stitch)<<"\n";

  if(!ok){
    // Localize the divergence. `pit_of` identifies the depression a cell belongs to
    // by that depression's pit cell (label-namespace-independent). raw (walk=false):
    // serial's wavefront leaf vs the resolved gLabel that feeds outlet discovery.
    // leaf (walk=true): the depression each cell ultimately resolves into.
    const auto pit_of = [&](const dh::DepressionHierarchy<float> &deps,
                            const rd::Array2D<dh_label_t> &lab, int x, int y, bool walk)->int{
      dh_label_t c = lab(x,y);
      // NOTE: lab is in the pre-collapse label namespace, so a cell whose leaf was compacted away
      // can carry a stale label >= deps.size(); bound the walk by range and depth so this diagnostic
      // never runs off the end (it can then only mis-report, never hang/UB).
      if(walk){ const float e=full(x,y);
        for(std::size_t g=0; c!=OCEAN && c<deps.size() && e>deps[c].out_elev && g<=deps.size(); g++)
          c=deps[c].parent; }
      if(c==OCEAN || c>=deps.size()) return -1;
      if(walk && deps[c].lchild!=dh::NO_VALUE) return -2;   // meta: pit cell is order-dependent
      return (int)deps[c].pit_cell;
    };
    int rawdiff=0, leafdiff=0;
    for(int y=0;y<full.height();y++) for(int x=0;x<full.width();x++){
      if(full.isNoData(x,y)) continue;
      if(pit_of(S,s_label,x,y,false)!=pit_of(G,gLabel,x,y,false)) rawdiff++;
      if(pit_of(S,s_label,x,y,true )!=pit_of(G,gLabel,x,y,true )) leafdiff++;
    }
    std::cerr<<"  raw-label diffs="<<rawdiff<<"  leaf-assignment diffs="<<leafdiff<<"\n";

    // For each cell whose resolved leaf differs, dump the depression chain in both
    // builds (label-namespace-independent: identify a node by its pit cell's x,y).
    const auto pxy = [&](const dh::DepressionHierarchy<float> &deps, dh_label_t c)->std::string{
      if(c==dh::NO_VALUE) return "-";
      if(c==OCEAN) return "OCEAN";
      if(c>=deps.size()) return "stale";                 // stale pre-collapse label (see note above)
      if(deps[c].pit_cell==dh::NO_VALUE) return "meta(no-pit)";
      int px,py; full.iToxy(deps[c].pit_cell,px,py); std::ostringstream os;
      os<<"pit("<<px<<","<<py<<")"; return os.str();
    };
    const auto dump_chain = [&](const char *tag, const dh::DepressionHierarchy<float> &deps,
                                const rd::Array2D<dh_label_t> &lab, int x, int y){
      std::cerr<<"    "<<tag<<" cell("<<x<<","<<y<<") elev="<<full(x,y)<<":\n";
      dh_label_t c = lab(x,y);
      for(int depth=0; c!=OCEAN && c<deps.size() && depth<12; depth++){
        const auto &d = deps[c];
        std::cerr<<"      "<<pxy(deps,c)<<" pit_elev="<<d.pit_elev<<" out_elev="<<d.out_elev
                 <<" lchild="<<pxy(deps,d.lchild)<<" rchild="<<pxy(deps,d.rchild)
                 <<" #ocean_linked="<<d.ocean_linked.size()<<"\n";
        c = d.parent;
      }
    };
    int shown=0;
    for(int y=0;y<full.height() && shown<3;y++) for(int x=0;x<full.width() && shown<3;x++){
      if(full.isNoData(x,y)) continue;
      if(pit_of(S,s_label,x,y,true)!=pit_of(G,gLabel,x,y,true)){
        std::cerr<<"  --- leaf-diff cell #"<<shown<<" ---\n";
        dump_chain("SERIAL", S, s_label, x, y);
        dump_chain("STITCH", G, gLabel, x, y);
        // Full field set of the stitch node the cell lands in, and who references it.
        const dh_label_t c = gLabel(x,y);
        if(c<G.size()){
          const auto &d = G[c];
          std::cerr<<"    STITCH node "<<c<<" "<<pxy(G,c)<<": parent="<<pxy(G,d.parent)
                   <<" odep="<<pxy(G,d.odep)<<" geolink="<<pxy(G,d.geolink)
                   <<" ocean_parent="<<d.ocean_parent<<" cc="<<d.cell_count<<" dv="<<d.dep_vol<<" ol={";
          for(auto o: d.ocean_linked) std::cerr<<pxy(G,o)<<" ";
          std::cerr<<"}\n";
        } else std::cerr<<"    STITCH node "<<c<<" (stale pre-collapse label)\n";
        shown++;
      }
    }

    std::cerr<<"  serial: "<<sig_serial.substr(0,300)<<"\n";
    std::cerr<<"  stitch: "<<sig_stitch.substr(0,300)<<"\n";
  }
  return (ok && flowdir_ok) ? 0 : 1;   // bit-identity = tree (canonical) AND per-cell flowdirs
}
