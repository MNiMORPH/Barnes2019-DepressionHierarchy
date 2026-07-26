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
#include <limits>
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

// §3.2 collapse pass -- contract seam-split artifacts into a serial-identical tree.
//
// A tiled flood can turn a cell (or flat) whose TRUE drainage exits across a seam
// into a spurious degenerate depression: a zero-height leaf (pit_elev==out_elev)
// whose pit sits on a tile edge with a cross-tile neighbour AT OR BELOW its elevation
// -- the escape the tile could not see (a lower neighbour => a monotonic slope cut by
// the seam; an equal neighbour => a flat straddling the seam, §3.2's deferred case).
// Both conditions are load-bearing:
//   * pit_elev==out_elev alone is NOT sufficient -- a serial flood CAN produce a
//     legitimate zero-height depression (a flat with an equal-elevation exit to a
//     neighbouring basin), which must be kept. (Observed: seed9 beta2.1.)
//   * the cross-tile "<=pit" escape is the seam-locality that distinguishes the
//     artifact (§3.2 "pit on a tile edge with a strictly-lower neighbour across",
//     here relaxed to <= to also fold the equal-elevation flat).
// A real depression's pit is a strict local minimum, so it never has a lower/equal
// cross-tile neighbour -- the test has no false positives on genuine basins.
//
// Such artifacts appear in two forms, contracted differently (Pass A below):
//   * ocean-linked splice -- the real basin lives in one tile as a node P; the seam
//     manufactured an extra degenerate leaf ocean_linked into P. Drop the artifact and
//     reattach its ocean_linked children to P. It holds no volume and its lone high
//     cell is above P's outlet, so P's cell_count/dep_vol already match serial.
//   * meta dissolve -- the basin's PIT straddles the seam, so it belongs to neither
//     tile; the stitch rebuilds it as a meta over the two tile-half artifact leaves.
//     That meta already carries the whole-basin aggregates, so dissolve it into one
//     leaf and drop both halves. (In this column-split harness a pit straddles exactly
//     one vertical seam -> two halves; a 2-D distributed build could split a corner
//     pit N ways, which would need the recursive subtree-dissolve generalisation.)
// Then labels are compacted. Grid-locality (pit cell, tile of a column) uses the DEM
// and tile bounds -- the 1-cell perimeter strips a distributed build already
// exchanges. Returns the number of artifacts contracted. O(#depressions) + O(boundary).
static int CollapseSeamArtifacts(dh::DepressionHierarchy<float> &G,
                                 const rd::Array2D<float> &full,
                                 const std::vector<int> &bounds){
  const dh_label_t N = G.size();
  std::vector<char> dead(N, 0);
  int contracted = 0, binary_skipped = 0;

  const auto tile_of = [&](int x){ int t=0; while(t+1<(int)bounds.size() && x>=bounds[t+1]) t++; return t; };

  // is_seam_artifact(i): is leaf i one the tiling manufactured? -- degenerate
  // (pit_elev==out_elev) with a cross-tile D8 neighbour at or below its pit (the
  // escape the local flood could not see). See the criterion above.
  const auto is_seam_artifact = [&](dh_label_t i)->bool {
    const auto &d = G[i];
    if(d.lchild!=dh::NO_VALUE || d.rchild!=dh::NO_VALUE) return false;  // leaves only
    if(d.pit_cell==dh::NO_VALUE) return false;                         // real pit, not a meta
    if(d.pit_elev!=d.out_elev) return false;                          // zero-height
    int px,py; full.iToxy(d.pit_cell, px, py);
    const float pe = d.pit_elev;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
      if(!dx && !dy) continue;
      const int nx=px+dx, ny=py+dy;
      if(!full.inGrid(nx,ny) || full.isNoData(nx,ny)) continue;
      if(tile_of(nx)!=tile_of(px) && full(nx,ny)<=pe) return true;
    }
    return false;
  };

  // Pass A -- mark artifacts for contraction. A seam-cut basin appears in one of two
  // forms depending on where its pit sits relative to the seam:
  //   * ocean-linked splice -- the real basin lives in one tile as a separate node P;
  //     the seam manufactured an extra degenerate leaf ocean_linked into P. Drop the
  //     artifact; its own ocean_linked children reattach to P.
  //   * meta dissolve -- the basin's PIT straddles the seam, so it has no home in
  //     either tile; the stitch rebuilds it as a meta over the two tile-half artifact
  //     leaves. That meta already carries the whole-basin aggregates (out_elev/
  //     cell_count/dep_vol == serial's), so dissolve it into one leaf and drop both
  //     halves. (Volume is conserved either way -- the refinement is exact.)
  for(dh_label_t i=1;i<N;i++){                        // skip ocean (node 0)
    if(dead[i]) continue;
    if(!is_seam_artifact(i)) continue;

    const dh_label_t P = G[i].parent;
    const auto &pol = G[P].ocean_linked;
    if(std::find(pol.begin(), pol.end(), i)!=pol.end()){
      dead[i] = 1;                                     // ocean-linked splice
      contracted++;
      continue;
    }

    // Binary child of meta P: dissolve P iff BOTH its children are seam artifacts
    // (the two halves of a pit-straddles-seam basin).
    const dh_label_t a = G[P].lchild, b = G[P].rchild;
    if(a!=dh::NO_VALUE && b!=dh::NO_VALUE && is_seam_artifact(a) && is_seam_artifact(b)){
      // The meta becomes the basin's single leaf. Its pit is the flood's seed cell --
      // the higher-index of the tied-lowest halves (the flood pops highest-index
      // first) -- at the shared floor elevation; its aggregates are already serial's.
      const dh_label_t keep = (G[a].pit_cell>=G[b].pit_cell) ? a : b;
      G[P].pit_cell = G[keep].pit_cell;
      G[P].pit_elev = std::min(G[a].pit_elev, G[b].pit_elev);
      G[P].lchild   = dh::NO_VALUE;
      G[P].rchild   = dh::NO_VALUE;
      dead[a] = 1; dead[b] = 1;
      contracted += 2;
      continue;
    }

    std::cerr<<"collapse: seam artifact "<<i<<" under meta "<<P
             <<" is not a two-halves dissolve (unhandled form); skipped\n";
    binary_skipped++;
  }
  (void)binary_skipped;
  if(contracted==0) return 0;

  // resolve(x) -- follow the parent chain up until a LIVE node. Artifacts are always
  // ocean_linked to their parent, so the parent chain is the ocean_linked spine; a
  // chain of stacked artifacts (a flat split by several seams) collapses to the first
  // live container. The ocean (node 0, never dead) terminates every chain.
  const auto resolve = [&](dh_label_t x)->dh_label_t {
    for(unsigned g=0; dead[x] && g<=N; g++) x = G[x].parent;
    return x;
  };

  // Compact: drop dead nodes, remap the survivors densely.
  std::vector<dh_label_t> perm(N, dh::NO_VALUE);
  dh_label_t next = 0;
  for(dh_label_t i=0;i<N;i++) if(!dead[i]) perm[i] = next++;
  const auto mv = [&](dh_label_t x)->dh_label_t {   // dead references redirect to the live container
    return x==dh::NO_VALUE ? dh::NO_VALUE : perm[resolve(x)];
  };

  dh::DepressionHierarchy<float> H(next);
  for(dh_label_t i=0;i<N;i++){
    if(dead[i]) continue;
    auto d = G[i];                                   // copy attributes
    d.parent    = mv(d.parent);
    d.odep      = mv(d.odep);
    d.geolink   = mv(d.geolink);
    d.lchild    = mv(d.lchild);                      // binary children are never artifacts
    d.rchild    = mv(d.rchild);
    d.dep_label = perm[i];
    d.ocean_linked.clear();                          // rebuilt below from the contracted edges
    H[perm[i]] = std::move(d);
  }
  // Rebuild the ocean_linked forest by contracting every original edge: a live child
  // re-homes to the first live node above its (possibly dead) parent; edges into dead
  // children vanish (that child's own children reconnect when its edges are processed).
  for(dh_label_t P=0;P<N;P++)
    for(const dh_label_t child : G[P].ocean_linked)
      if(!dead[child])
        H[perm[resolve(P)]].ocean_linked.push_back(perm[child]);

  G = std::move(H);
  return contracted;
}

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

  // ---- per-tile: pre-label seam-crossing edges BOUNDARY, then PhaseAB ----
  // gFlow is the one authoritative drainage map: the flood's own flowdirs, assembled
  // globally. Everything downhill (conduit resolution) follows it, so drainage is
  // internally consistent with the flood -- including flats, where the flood's
  // deterministic order (radix_heap bucket sort) fixes the direction.
  rd::Array2D<int8_t> gFlow(full.width(), full.height(), rd::NO_FLOW);
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
    for(int y=0;y<tile.dem.height();y++)                 // assemble the global flowdir map
      for(int x=0;x<tile.dem.width();x++)
        gFlow(tile.x0 + x, y) = fd(x, y);
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
  // into, by FOLLOWING THE FLOOD'S FLOWDIRS (gFlow) -- the same drainage map as the
  // flood, so it is flat-consistent by construction -- all the way to a real terminal
  // (a depression or the ocean). A seam-edge BOUNDARY seed has NO_FLOW (its downstream
  // is across the seam), so there we cross via the lowest cross-seam neighbour. Order-
  // independent: resolutions are computed from the original labels, then applied. ----
  {
    const auto terminal = [&](int x, int y)->dh_label_t {
      for(unsigned int guard=0; guard<=full.size(); guard++){
        const auto l = gLabel(x,y);
        if(l!=BOUNDARY) return l;                          // reached a depression or ocean
        const int fd = gFlow(x,y);
        if(fd!=rd::NO_FLOW){ x += rd::d8x[fd]; y += rd::d8y[fd]; }
        else {                                             // seam seed: cross the seam
          int lx,ly; bool to_ocean;
          if(!drain(x,y,lx,ly,to_ocean) || to_ocean) return OCEAN;
          x=lx; y=ly;
        }
      }
      return OCEAN;                                        // safety net (no cycle expected)
    };
    std::vector<std::pair<int,int>> bcells;
    for(int y=0;y<full.height();y++)
      for(int x=0;x<full.width();x++)
        if(gLabel(x,y)==BOUNDARY) bcells.emplace_back(x,y);
    std::vector<dh_label_t> resolved(bcells.size());
    for(size_t i=0;i<bcells.size();i++) resolved[i] = terminal(bcells[i].first, bcells[i].second);
    for(size_t i=0;i<bcells.size();i++) gLabel(bcells[i].first, bcells[i].second) = resolved[i];
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

  // ---- §3.2 collapse pass: contract seam-split artifacts to a serial-identical tree ----
  const int n_collapsed = CollapseSeamArtifacts(G, full, bounds);
  if(n_collapsed) std::cerr<<"collapse: contracted "<<n_collapsed<<" seam artifact(s)\n";

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

    // For each cell whose resolved leaf differs, dump the depression chain in both
    // builds (label-namespace-independent: identify a node by its pit cell's x,y).
    const auto pxy = [&](const dh::DepressionHierarchy<float> &deps, dh_label_t c)->std::string{
      if(c==dh::NO_VALUE) return "-";
      if(c==OCEAN) return "OCEAN";
      if(deps[c].pit_cell==dh::NO_VALUE) return "meta(no-pit)";
      int px,py; full.iToxy(deps[c].pit_cell,px,py); std::ostringstream os;
      os<<"pit("<<px<<","<<py<<")"; return os.str();
    };
    const auto dump_chain = [&](const char *tag, const dh::DepressionHierarchy<float> &deps,
                                const rd::Array2D<dh_label_t> &lab, int x, int y){
      std::cerr<<"    "<<tag<<" cell("<<x<<","<<y<<") elev="<<full(x,y)<<":\n";
      dh_label_t c = lab(x,y);
      for(int depth=0; c!=OCEAN && depth<12; depth++){
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
        const auto &d = G[c];
        std::cerr<<"    STITCH node "<<c<<" "<<pxy(G,c)<<": parent="<<pxy(G,d.parent)
                 <<" odep="<<pxy(G,d.odep)<<" geolink="<<pxy(G,d.geolink)
                 <<" ocean_parent="<<d.ocean_parent<<" cc="<<d.cell_count<<" dv="<<d.dep_vol<<" ol={";
        for(auto o: d.ocean_linked) std::cerr<<pxy(G,o)<<" ";
        std::cerr<<"}\n";
        shown++;
      }
    }

    std::cerr<<"  serial: "<<sig_serial.substr(0,300)<<"\n";
    std::cerr<<"  stitch: "<<sig_stitch.substr(0,300)<<"\n";
  }
  return ok ? 0 : 1;
}
