// dephier_stitch -- in-process tiled DepressionHierarchy build, diffed against serial.
//
// Stitch with BOUNDARY seeding (PARALLEL_DEPHIER_PLAN.md §3, §7.1):
//
//   per tile:  extract sub-DEM -> label ocean; pre-label seam-edge cells whose
//              steepest descent CROSSES the seam as BOUNDARY (so they seed as
//              provisional exterior and never become spurious pits) -> PhaseAB
//   global:    remap tiles into one namespace (shared ocean, offset the rest);
//              resolve every BOUNDARY cell to the real depression it drains into
//              (steepest descent, lowest first); build the outlet set as per-tile
//              intra-tile outlets + Barnes' HandleEdge across the seam perimeters;
//              run one global PhaseCD.
//   verify:    canonical signature must equal the serial build's.
//
// The outlet step is in the distributable shape: each tile derives outlets from its
// own resolved labels, and only the 1-cell perimeter strips cross the seam (Barnes'
// HandleEdge). In this in-process harness that is organizationally equivalent to a
// single global pass; in a distributed build it is what keeps the per-rank footprint
// O(N/P) + O(boundary). BOUNDARY conduit resolution still uses the full grid here
// (the distributed version needs a cross-tile boundary-graph pass -- future work).

#include "dh_canonical.hpp"

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
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
  dh::DepressionHierarchy<elev_t> deps;
  std::vector<dh::Outlet<elev_t>> outlets;
  int        x0 = 0;
  dh_label_t offset = 0;
  dh_label_t g(dh_label_t local) const {
    if(local==OCEAN || local==BOUNDARY) return local;   // sentinels are global as-is
    return offset + (local - 1);
  }
};

int main(int argc, char **argv){
  if(argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level> <Split Cols (comma-sep)>\n";
    return -1;
  }
  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);

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
  // Lowest LAND neighbour on the full grid, and whether the cell touches ocean.
  // An ocean neighbour is the sea (effectively -inf), so a cell touching ocean
  // drains to the ocean regardless of its land neighbours.
  const auto drain = [&](int x,int y,int &lx,int &ly,bool &to_ocean)->bool{
    float best = full(x,y); bool found=false; to_ocean=false;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
      if(!dx && !dy) continue;
      int nx=x+dx, ny=y+dy;
      if(!full.inGrid(nx,ny)) continue;
      if(is_ocean(nx,ny)){ to_ocean=true; continue; }
      if(full(nx,ny) < best){ best=full(nx,ny); lx=nx; ly=ny; found=true; }
    }
    return found;
  };

  // ---- per-tile: pre-label seam-crossing edges BOUNDARY, then PhaseAB ----
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

    rd::Array2D<int8_t> fd(tile.dem.width(), tile.dem.height(), rd::NO_FLOW);
    dh::GetDepressionHierarchyPhaseAB<float,rd::Topology::D8>(tile.dem, tile.label, fd, tile.deps, tile.outlets);
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

  // ---- conduit resolution: every BOUNDARY cell adopts the depression it drains
  // into. Resolve lowest-first, so each cell's steepest-descent neighbour (which is
  // strictly lower) is already a real label. ----
  {
    std::vector<std::pair<float,int>> bcells;   // (elev, flat index)
    for(int y=0;y<full.height();y++)
      for(int x=0;x<full.width();x++)
        if(gLabel(x,y)==BOUNDARY)
          bcells.emplace_back(full(x,y), full.xyToI(x,y));
    std::sort(bcells.begin(), bcells.end());
    for(auto &bc : bcells){
      int x,y; full.iToxy(bc.second, x, y);
      int lx,ly; bool to_ocean;
      const bool found = drain(x,y,lx,ly,to_ocean);
      if(to_ocean || !found) gLabel(x,y) = OCEAN;          // drains to the sea
      else                   gLabel(x,y) = gLabel(lx,ly);  // drains into a depression
    }
  }

  // ---- global outlet set, in the distributable shape (Barnes' join): each tile
  // derives its own outlets from its resolved labels (intra-tile adjacencies), and
  // only the perimeter strips cross the seam via HandleEdge. Both feed one database
  // keyed on the depression pair, keeping the lowest max-of-pair -- the serial PhaseB
  // rule. Intra-tile + cross-seam together cover every adjacency, so the result is
  // identical to a single global pass. ----
  std::vector<dh::Outlet<float>> outlets;
  {
    std::map<std::pair<dh_label_t,dh_label_t>, std::pair<float,dh::flat_c_idx>> db;
    const auto record = [&](dh_label_t la, dh_label_t lb,
                            float ea, dh::flat_c_idx ca, float eb, dh::flat_c_idx cb){
      if(la==lb) return;                          // same depression: not an outlet
      float oelev; dh::flat_c_idx ocell;          // outlet = the higher of the pair
      if(ea>=eb){ oelev=ea; ocell=ca; } else { oelev=eb; ocell=cb; }
      const auto key = std::minmax(la, lb);
      const auto it = db.find({key.first,key.second});
      if(it==db.end() || oelev < it->second.first)
        db[{key.first,key.second}] = {oelev, ocell};
    };

    // Per-tile: adjacencies whose neighbour lies in the SAME tile.
    for(int y=0;y<full.height();y++)
      for(int x=0;x<full.width();x++){
        if(full.isNoData(x,y)) continue;
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
          if(!dx && !dy) continue;
          int nx=x+dx, ny=y+dy;
          if(!full.inGrid(nx,ny) || full.isNoData(nx,ny)) continue;
          if(tile_of(nx)!=tile_of(x)) continue;   // cross-seam: handled below
          record(gLabel(x,y), gLabel(nx,ny), full(x,y), full.xyToI(x,y),
                                             full(nx,ny), full.xyToI(nx,ny));
        }
      }

    // HandleEdge: match the two resolved perimeter strips at each seam (D8), exactly
    // as Barnes' parallel priority-flood pairs edge cell i with neighbours i-1,i,i+1.
    for(size_t b=1;b+1<bounds.size();b++){
      const int cA=bounds[b]-1, cB=bounds[b];     // the two columns straddling the seam
      for(int y=0;y<full.height();y++){
        if(full.isNoData(cA,y)) continue;
        for(int ny=y-1;ny<=y+1;ny++){
          if(ny<0 || ny>=full.height() || full.isNoData(cB,ny)) continue;
          record(gLabel(cA,y), gLabel(cB,ny), full(cA,y), full.xyToI(cA,y),
                                              full(cB,ny), full.xyToI(cB,ny));
        }
      }
    }

    for(auto &kv : db){
      dh::Outlet<float> o;
      o.depa = kv.first.first; o.depb = kv.first.second;
      o.out_elev = kv.second.first; o.out_cell = kv.second.second;
      outlets.push_back(o);
    }
  }

  // ---- one global PhaseCD ----
  dh::GetDepressionHierarchyPhaseCD<float>(G, outlets, full, gLabel);

  // ---- serial ground truth ----
  auto s_label = ocean_labels(full, ocean_level);
  rd::Array2D<int8_t> s_fd(full.width(), full.height(), rd::NO_FLOW);
  auto S = dh::GetDepressionHierarchy<float,rd::Topology::D8>(full, s_label, s_fd);

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

  if(!ok){
    // Localize the divergence. `pit_of` identifies the depression a cell belongs to
    // by that depression's pit cell (label-namespace-independent). raw (walk=false):
    // serial's wavefront leaf vs the resolved gLabel that feeds outlet discovery.
    // leaf (walk=true): the depression each cell ultimately resolves into.
    const auto pit_of = [&](const dh::DepressionHierarchy<float> &deps,
                            const rd::Array2D<dh_label_t> &lab, int x, int y, bool walk)->int{
      dh_label_t c = lab(x,y);
      if(walk){ const float e=full(x,y); while(c!=OCEAN && e>deps[c].out_elev) c=deps[c].parent; }
      if(c==OCEAN) return -1;
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
    std::cerr<<"  serial: "<<sig_serial.substr(0,300)<<"\n";
    std::cerr<<"  stitch: "<<sig_stitch.substr(0,300)<<"\n";
  }
  return ok ? 0 : 1;
}
