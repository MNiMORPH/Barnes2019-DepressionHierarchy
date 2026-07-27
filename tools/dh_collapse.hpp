// dh_collapse -- the §3.2 collapse pass, shared by the in-process stitch
// (tools/dephier_stitch.cpp) and the distributed harness (tools/dephier_mpi.cpp) so the
// seam-artifact contraction is defined once. Extracted verbatim from the stitch; no behaviour
// change. See the comment on CollapseSeamArtifacts below.
#pragma once

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;

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
inline int CollapseSeamArtifacts(dh::DepressionHierarchy<float> &G,
                                 const rd::Array2D<float> &full,
                                 const std::vector<int> &bounds){
  using dh::dh_label_t;
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
