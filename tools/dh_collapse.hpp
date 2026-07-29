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
// Such artifacts appear in two forms, contracted differently (Passes A and B below):
//   * ocean-linked splice (Pass A) -- the real basin lives in one tile as a node P; the seam
//     manufactured an extra degenerate leaf, seam-adjacent (is_seam_artifact) and ocean_linked
//     into P. Drop the artifact and reattach its ocean_linked children to P. It holds no volume
//     and its lone high cell is above P's outlet, so P's cell_count/dep_vol already match serial.
//   * meta dissolve (Pass B) -- the basin's FLOOR straddles the seam, so the flat belongs to
//     neither tile; the stitch rebuilds it as a meta over the two tile-half leaves, which are
//     degenerate leaves whose pits fall in DIFFERENT tiles. (This straddle test, not the pit's
//     seam-adjacency, is what identifies the halves -- a flat's fragment pit can land in the tile
//     interior.) That meta already carries the whole-basin aggregates, so dissolve it into one
//     leaf and drop both halves. (In this column-split harness a flat straddles exactly one
//     vertical seam -> two halves; a 2-D distributed build could split a corner flat N ways,
//     which would need the recursive subtree-dissolve generalisation.)
// Then labels are compacted. Grid-locality (pit cell, tile of a column) uses the DEM
// and tile bounds -- the 1-cell perimeter strips a distributed build already
// exchanges. Returns the number of artifacts contracted. O(#depressions) + O(boundary).
inline int CollapseSeamArtifacts(dh::DepressionHierarchy<float> &G,
                                 const rd::Array2D<float> &full,
                                 const std::vector<int> &bounds){
  using dh::dh_label_t;
  const dh_label_t N = G.size();
  std::vector<char> dead(N, 0);
  int contracted = 0;
  int meta_over_halves = 0;   // Pass B / B2 firings: the seam-DEPENDENT dissolutions we want to RETIRE
                              // (a seam-split flat rebuilt as a meta-over-halves, dissolved with a
                              // seam-chosen pit). Early flat unification (split-invariant-flats mode)
                              // should drive this to zero; a nonzero count is warned below.

  const auto tile_of = [&](int x){ int t=0; while(t+1<(int)bounds.size() && x>=bounds[t+1]) t++; return t; };

  // is_seam_artifact(i): is leaf i one the tiling manufactured? -- degenerate
  // (pit_elev==out_elev) with a cross-tile D8 neighbour at or below its pit (the escape the local
  // flood could not see). Used ONLY by Pass A (ocean-linked splice); Pass B uses a straddle test.
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
      if(!full.inGrid(nx,ny) || tile_of(nx)==tile_of(px)) continue;    // only cross-seam escapes
      // A cross-tile NoData cell is OCEAN (base level, below any land pit) -- the escape the tile could
      // not see; a cross-tile land cell at/below the pit is likewise an escape. Reading NoData as OCEAN
      // here matches ocean_labels and the record() outlet scans (the collapse has no label grid). Without
      // it, a coastal degenerate pit whose only lower neighbour is NoData-ocean across the seam is missed
      // -- the split-dependent extra zero-volume leaf (e.g. Corsica pit(78,63) at split 78).
      if(full.isNoData(nx,ny) || full(nx,ny)<=pe) return true;
    }
    return false;
  };

  // A real zero-height leaf, regardless of where its pit sits relative to a seam. (The
  // meta-dissolve reads seam-membership from the two halves straddling the seam, not from the
  // pit's seam-adjacency, so it does NOT use is_seam_artifact -- see Pass B.)
  const auto is_degenerate_leaf = [&](dh_label_t i)->bool {
    const auto &d = G[i];
    return d.lchild==dh::NO_VALUE && d.rchild==dh::NO_VALUE   // leaf
        && d.pit_cell!=dh::NO_VALUE                           // real pit, not the ocean
        && d.pit_elev==d.out_elev;                            // zero-height
  };
  const auto pit_tile = [&](dh_label_t i){ int x,y; full.iToxy(G[i].pit_cell,x,y); return tile_of(x); };

  // Pass A -- ocean-linked splice. The real basin lives in one tile as a node P; the seam
  // manufactured an extra degenerate leaf, seam-adjacent (is_seam_artifact) and ocean_linked into
  // P. Drop the artifact; its own ocean_linked children reattach to the live container via resolve()
  // during compaction. (This is the form where the basin's pit did NOT straddle the seam.)
  for(dh_label_t i=1;i<N;i++){                        // skip ocean (node 0)
    if(dead[i]) continue;
    if(!is_seam_artifact(i)) continue;
    const dh_label_t P = G[i].parent;
    const auto &pol = G[P].ocean_linked;
    if(std::find(pol.begin(), pol.end(), i)!=pol.end()){
      dead[i] = 1;
      contracted++;
    }
  }

  // Pass B -- meta dissolve. A meta whose two children are both degenerate leaves whose pits lie in
  // DIFFERENT tiles is a single serial flat basin the seam split into halves: the basin's floor
  // straddles the seam, so it has no home in either tile and the stitch rebuilt it as a meta over the
  // two tile-half leaves. Seam-membership is read from the halves straddling the seam, NOT the pit's
  // seam-adjacency -- so this also catches a flat whose fragment pit landed in the tile interior (the
  // pit-only is_seam_artifact test missed exactly that; see ENH-2). The meta already carries the
  // whole-basin aggregates (out_elev/cell_count/dep_vol == serial's), so dissolve it into one leaf and
  // drop both halves; volume is conserved, the refinement is exact. Ascending order dissolves nested
  // metas (a flat cut by several seams) bottom-up. Two zero-volume degenerate leaves in different
  // tiles under one meta cannot arise without a seam, so there are no false positives on genuine basins.
  for(dh_label_t P=1;P<N;P++){                        // skip ocean (node 0)
    if(dead[P]) continue;
    const dh_label_t a = G[P].lchild, b = G[P].rchild;
    if(a==dh::NO_VALUE || b==dh::NO_VALUE) continue;  // P is not a meta
    if(dead[a] || dead[b]) continue;
    if(!is_degenerate_leaf(a) || !is_degenerate_leaf(b)) continue;
    if(pit_tile(a)==pit_tile(b)) continue;            // the two halves must straddle a seam
    // The meta becomes the basin's single leaf. Its pit is the flood's seed cell -- the higher-index
    // of the tied-lowest halves (the flood pops highest-index first) -- at the shared floor elevation.
    const dh_label_t keep = (G[a].pit_cell>=G[b].pit_cell) ? a : b;
    G[P].pit_cell = G[keep].pit_cell;
    G[P].pit_elev = std::min(G[a].pit_elev, G[b].pit_elev);
    G[P].lchild   = dh::NO_VALUE;
    G[P].rchild   = dh::NO_VALUE;
    dead[a] = 1; dead[b] = 1;
    contracted += 2; meta_over_halves++;
  }

  // Pass B2 -- rim-fragment dissolve. Pass B handles a basin whose FLOOR straddles a seam (both halves
  // degenerate). This is the complementary case: a basin whose floor lies wholly in ONE tile but whose RIM
  // the seam cut. The tiling manufactures the across-seam rim strip as a spurious degenerate "depression"
  // (a seam artifact) that merges back at the rim elevation, so the stitch rebuilds the single serial basin
  // as a META over {the real basin leaf, the degenerate rim fragment}. Serial keeps it as ONE leaf (the real
  // basin, outletting past the rim). Detect: a meta with EXACTLY ONE is_seam_artifact child (the rim
  // fragment -- degenerate with a cross-tile escape into its sibling) and one NON-degenerate real leaf.
  // Dissolve into the real child's pit, keeping the meta's out_elev/cell_count/dep_vol -- which already equal
  // serial's whole basin (volume conserved, refinement exact). The is_seam_artifact guard (strict-local-min
  // pits never have a lower/equal cross-tile neighbour) is what keeps this off genuine basins; the real side
  // being non-degenerate excludes the Pass B floor-straddle. (Observed: testdem8 split 3 -- the 4-rim of the
  // 2-basin cut at the seam; serial 1 leaf (2,9,36,212) vs stitch meta over (2,4,30,32)+(4,4,6,0).)
  for(dh_label_t P=1;P<N;P++){                        // skip ocean (node 0)
    if(dead[P]) continue;
    const dh_label_t a = G[P].lchild, b = G[P].rchild;
    if(a==dh::NO_VALUE || b==dh::NO_VALUE) continue;  // P is not a meta
    if(dead[a] || dead[b]) continue;
    const bool a_art = is_seam_artifact(a), b_art = is_seam_artifact(b);
    if(a_art == b_art) continue;                       // need EXACTLY one seam-artifact child
    const dh_label_t art  = a_art ? a : b;             // the spurious degenerate rim fragment
    const dh_label_t real = a_art ? b : a;             // the genuine basin the rim belongs to
    if(G[real].lchild!=dh::NO_VALUE || G[real].rchild!=dh::NO_VALUE) continue;  // real side must be a leaf
    if(is_degenerate_leaf(real)) continue;             // ...and a genuine (non-degenerate) basin
    G[P].pit_cell = G[real].pit_cell;                  // the meta becomes the real basin's single leaf,
    G[P].pit_elev = G[real].pit_elev;                  // keeping its own out_elev/cell_count/dep_vol
    G[P].lchild   = dh::NO_VALUE;
    G[P].rchild   = dh::NO_VALUE;
    dead[art] = 1; dead[real] = 1;
    contracted += 2; meta_over_halves++;
  }

  // Retirement warning: a meta-over-halves dissolution is seam-DEPENDENT (its pit is a seam-chosen
  // representative, and whether it even fires depends on how the seam cut the flat -- so it does a
  // different amount of work per tiling). It is a bridge we want to remove by unifying seam-split flats
  // EARLY (before the meta forms); when that mode fully catches them this count is 0. See
  // SPLIT_INVARIANT_FLATS_PLAN.md.
  if(meta_over_halves>0)
    std::cerr<<"WARNING: collapse dissolved "<<meta_over_halves<<" meta-over-halves seam artifact(s) "
             <<"(seam-dependent; should retire via early split-invariant flat unification)\n";

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
