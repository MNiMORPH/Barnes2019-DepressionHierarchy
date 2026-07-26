// dephier_stitch -- in-process tiled DepressionHierarchy build, diffed against serial.
//
// The first realization of the stitch (PARALLEL_DEPHIER_PLAN.md §7.1) on the benign
// case: a domain with an ocean ring, split into 1xN column tiles, so every tile
// touches ocean and no BOUNDARY seeding is needed yet (that arrives with interior
// tiles). It reuses the exact serial PhaseAB / PhaseCD (no dephier.hpp changes):
//
//   per tile:  extract sub-DEM -> label ocean -> PhaseAB  => {depressions, outlets}
//   global:    remap tiles into one namespace (shared ocean=0, offset the rest)
//              match perimeter strips at each seam (Barnes' HandleEdge rule) -> cross outlets
//              assemble global {depressions, outlets, label grid}
//              run ONE global PhaseCD => global tree
//   verify:    canonical signature (dh_canonical) must equal the serial build's.
//
// Because the benign fixtures are tie-free, the match is bit-identity.

#include "dh_canonical.hpp"

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

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

// Columns [x0,x1) of `full` as a standalone tile, preserving the NoData value.
static rd::Array2D<float> extract_cols(const rd::Array2D<float> &full, int x0, int x1){
  rd::Array2D<float> t(x1 - x0, full.height());
  t.setNoData(full.noData());
  for(int y=0;y<full.height();y++)
    for(int x=x0;x<x1;x++)
      t(x - x0, y) = full(x, y);
  return t;
}

// Label ocean (NoData or == ocean_level) as OCEAN, everything else NO_DEP.
static rd::Array2D<dh_label_t> ocean_labels(const rd::Array2D<float> &dem, float ocean_level){
  rd::Array2D<dh_label_t> label(dem.width(), dem.height(), NO_DEP);
  for(unsigned int i=0;i<dem.size();i++)
    if(dem.isNoData(i) || dem(i)==ocean_level)
      label(i) = OCEAN;
  return label;
}

// One tile's PhaseAB result plus its local->global label map.
template<class elev_t>
struct Tile {
  rd::Array2D<elev_t>          dem;
  rd::Array2D<dh_label_t>      label;    // leaf labels after PhaseAB (local namespace)
  dh::DepressionHierarchy<elev_t> deps;  // [ocean, leaves...] (local)
  std::vector<dh::Outlet<elev_t>> outlets;
  int x0 = 0;                            // global column of this tile's left edge
  dh_label_t offset = 0;                 // global(local k>0) = offset + (k-1)
  dh_label_t g(dh_label_t local) const { return local==OCEAN ? OCEAN : offset + (local - 1); }
};

int main(int argc, char **argv){
  if(argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level> <Split Cols (comma-sep)>\n";
    std::cout<<"  e.g. "<<argv[0]<<" int_two.dem -9999 20     (two tiles: [0,20) and [20,W))\n";
    return -1;
  }
  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);

  // Parse split columns into tile boundaries.
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

  // ---- per-tile PhaseAB ----
  std::vector<Tile<float>> tiles(ntiles);
  dh_label_t next_offset = 1;  // global label 0 is the shared ocean
  for(int t=0;t<ntiles;t++){
    auto &tile = tiles[t];
    tile.x0    = bounds[t];
    tile.dem   = extract_cols(full, bounds[t], bounds[t+1]);
    tile.label = ocean_labels(tile.dem, ocean_level);
    rd::Array2D<int8_t> fd(tile.dem.width(), tile.dem.height(), rd::NO_FLOW);
    dh::GetDepressionHierarchyPhaseAB<float,rd::Topology::D8>(tile.dem, tile.label, fd, tile.deps, tile.outlets);
    tile.offset = next_offset;
    next_offset += tile.deps.size() - 1;   // number of non-ocean depressions in this tile
  }
  const dh_label_t n_global = next_offset;  // total global depression count

  // ---- assemble global depressions (shared ocean at 0, leaves remapped) ----
  dh::DepressionHierarchy<float> G(n_global);
  G[0] = tiles[0].deps[0];        // ocean prototype (PhaseAB leaves it empty of links)
  G[0].dep_label = 0;
  for(auto &tile : tiles)
    for(dh_label_t k=1;k<tile.deps.size();k++){
      auto d = tile.deps[k];
      d.dep_label = tile.g(k);
      if(d.pit_cell!=dh::NO_VALUE){          // globalise the tile-local pit index
        int lx,ly; tile.dem.iToxy(d.pit_cell, lx, ly);
        d.pit_cell = full.xyToI(tile.x0 + lx, ly);
      }
      G[tile.g(k)] = d;
    }

  // ---- assemble global outlets: intra-tile (remapped) ----
  std::vector<dh::Outlet<float>> outlets;
  for(auto &tile : tiles)
    for(auto o : tile.outlets){
      o.depa = tile.g(o.depa);
      o.depb = tile.g(o.depb);
      outlets.push_back(o);
    }

  // ---- cross-tile outlets: match adjacent perimeter strips (Barnes HandleEdge) ----
  std::map<std::pair<dh_label_t,dh_label_t>, float> cross;
  for(int t=0;t+1<ntiles;t++){
    auto &A = tiles[t];      // left tile: its right column touches the seam
    auto &B = tiles[t+1];    // right tile: its left column touches the seam
    const auto ea = A.dem.rightColumn();
    const auto eb = B.dem.leftColumn();
    const auto la = A.label.rightColumn();
    const auto lb = B.label.leftColumn();
    const int H = ea.size();
    for(int y=0;y<H;y++){
      const dh_label_t ga = A.g(la[y]);
      for(int ny=y-1;ny<=y+1;ny++){       // D8 across the seam
        if(ny<0 || ny>=H) continue;
        const dh_label_t gb = B.g(lb[ny]);
        if(ga==gb) continue;              // same label (e.g. both ocean): no outlet
        const float elev = std::max(ea[y], eb[ny]);
        auto key = std::minmax(ga, gb);
        auto it = cross.find({key.first, key.second});
        if(it==cross.end() || elev < it->second)
          cross[{key.first, key.second}] = elev;
      }
    }
  }
  for(auto &kv : cross){
    dh::Outlet<float> o;
    o.depa = kv.first.first;
    o.depb = kv.first.second;
    o.out_elev = kv.second;
    o.out_cell = dh::NO_VALUE;
    outlets.push_back(o);
  }

  // ---- global label grid for the volume pass ----
  rd::Array2D<dh_label_t> gLabel(full.width(), full.height(), OCEAN);
  for(auto &tile : tiles)
    for(int y=0;y<tile.dem.height();y++)
      for(int x=0;x<tile.dem.width();x++)
        gLabel(tile.x0 + x, y) = tile.g(tile.label(x, y));

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

  std::cerr<<"tiles="<<ntiles
           <<"  serial nodes="<<iv_serial.n_nodes<<" (leaf "<<iv_serial.n_leaf<<", meta "<<iv_serial.n_meta<<")"
           <<"  stitch nodes="<<iv_stitch.n_nodes<<" (leaf "<<iv_stitch.n_leaf<<", meta "<<iv_stitch.n_meta<<")\n";
  std::cerr<<"serial total_dep_vol="<<iv_serial.total_dep_vol<<"  stitch total_dep_vol="<<iv_stitch.total_dep_vol<<"\n";

  // Diagnostic: how many cells land in a different physical depression (identified
  // by its global pit cell) between serial and stitch, and in which columns?
  {
    // Resolve a cell to the depression that actually CONTAINS it: walk up from the
    // wavefront leaf label while the cell is above that depression's outlet (this is
    // what CalculateMarginalVolumes does). Identify it by its global pit cell.
    // -1 = resolves to ocean; -2 = resolves to a meta-depression (whose pit_cell is
    // an order-dependent choice of one child, so not comparable); else the leaf's
    // global pit cell. This isolates GENUINE leaf-membership differences.
    const auto resolved_pit = [&](const dh::DepressionHierarchy<float> &deps,
                                  const rd::Array2D<dh_label_t> &lab, int x, int y) -> int {
      dh_label_t c = lab(x,y);
      const float e = full(x,y);
      while(c!=OCEAN && e > deps[c].out_elev) c = deps[c].parent;
      if(c==OCEAN) return -1;
      if(deps[c].lchild!=dh::NO_VALUE) return -2;   // meta
      return (int)deps[c].pit_cell;
    };
    std::map<int,int> diff_by_col;
    int ndiff = 0, shown = 0;
    for(int y=0;y<full.height();y++)
      for(int x=0;x<full.width();x++){
        const auto spit = resolved_pit(S, s_label, x, y);
        const auto gpit = resolved_pit(G, gLabel, x, y);
        if(spit!=gpit){
          ndiff++; diff_by_col[x]++;
          if(shown++ < 8)
            std::cerr<<"  DIFF cell("<<x<<","<<y<<") elev="<<full(x,y)<<" s_pit="<<spit<<" g_pit="<<gpit<<"\n";
        }
      }
    std::cerr<<"partition diff: "<<ndiff<<" leaf cells assigned to a different depression";
    if(ndiff){ std::cerr<<" at columns"; for(auto &c:diff_by_col) std::cerr<<" "<<c.first<<"(x"<<c.second<<")"; }
    std::cerr<<"\n";
  }

  const bool ok = (sig_stitch==sig_serial);
  std::cout<<(ok ? "STITCH-MATCH " : "STITCH-DIFFER ")<<in_name<<" splits="<<argv[3]<<"\n";
  if(!ok){
    std::cerr<<"  signatures differ (stitch "<<sig_stitch.size()<<" chars, serial "<<sig_serial.size()<<")\n";
    std::cerr<<"  serial: "<<sig_serial<<"\n";
    std::cerr<<"  stitch: "<<sig_stitch<<"\n";
  }
  return ok ? 0 : 1;
}
