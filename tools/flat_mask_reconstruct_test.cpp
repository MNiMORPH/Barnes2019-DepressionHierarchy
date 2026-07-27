// flat_mask_reconstruct_test -- documents and verifies, by execution, the equivalence that the
// distributed flat resolution (ENH-1, see ENHANCEMENTS.md) is built on.
//
// THE EQUIVALENCE (three formulations of Barnes-2014's flat_mask that all agree):
//
//   (A) richdem's resolve_flats_barnes: labels each flat with label_this (a flood-fill over
//       equal-elevation D8-neighbours), then two gradient BFSs that only cross cells of the
//       SAME LABEL, and flat_height[label] = max away-distance in that label.
//
//   (B) LABEL-FREE reconstruction: because label_this floods equal-elevation D8-neighbours,
//       "same flat label" is IDENTICAL to "D8-adjacent AND same elevation." So the label check
//       in each BFS can be replaced by an elevation check, and flat_height becomes a plain
//       max-relaxation of the away-distance over same-elevation adjacency. No labels anywhere:
//         away_dist    = multi-source min-BFS from HIGH edges over {NO_FLOW, same-elevation}
//         towards_dist = multi-source min-BFS from LOW  edges over {NO_FLOW, same-elevation}
//         flat_height  = max-relaxation of away_dist   over {same-elevation}
//         flat_mask    = flat_height - away_dist + 2*towards_dist   (then d8_masked_FlowDir)
//
//   (C) DISTRIBUTED (ENH-1): each of the three fields in (B) is an ORDER-INDEPENDENT relaxation,
//       so it distributes as local-compute -> 1-column seam exchange -> relax -> iterate to a
//       global "changed?" all-reduce, and converges to the exact same fixed point as (B)==(A).
//       O(boundary) per round, rounds ~ flat diameter -- the memory-for-latency trade.
//
// This test proves (A) == (B) bit-for-bit on a synthetic flat. (C) is validated separately, in
// the harness, against (B). If this test ever fails, the label<->elevation-adjacency equivalence
// that ENH-1 relies on has been broken and the distributed flat resolution is unsound.
#include "dh_flats.hpp"   // resolve_flat_flowdirs (reference) + resolve_flat_flowdirs_option2 (ENH-1)

#include <richdem/common/Array2D.hpp>
#include <richdem/common/grid_cell.hpp>
#include <richdem/flowmet/d8_flowdirs.hpp>
#include <richdem/flats/flat_resolution.hpp>

#include <deque>
#include <iostream>
#include <vector>

namespace rd = richdem;
using rd::NO_FLOW;

int main(){
  // Synthetic "bathtub": high rim (10), flat floor (5), one drain to the edge -- a flat with
  // many high edges (floor against rim) and a low edge (floor cell beside the drain).
  const int W=24, H=16;
  rd::Array2D<float> dem(W,H,5.0f);
  for(int x=0;x<W;x++){ dem(x,0)=10; dem(x,H-1)=10; }
  for(int y=0;y<H;y++){ dem(0,y)=10; dem(W-1,y)=10; }
  dem(12,H-1)=0;                        // outlet
  dem(12,H-2)=2;                        // small ramp so the floor drains to it

  // (A) richdem reference.
  rd::Array2D<int8_t> sd(W,H,NO_FLOW);
  rd::d8_flow_directions(dem, sd);
  rd::Array2D<int32_t> fm_ref, lab;
  rd::resolve_flats_barnes(dem, sd, fm_ref, lab);

  // (B) label-free reconstruction over same-elevation adjacency.
  // High/low edges exactly as find_flat_edges: a LOW edge has flow and touches a same-elevation
  // NO_FLOW cell; a HIGH edge is NO_FLOW and touches a higher cell.
  std::deque<rd::GridCell> high_edges, low_edges;
  for(int x=0;x<W;x++) for(int y=0;y<H;y++)
    for(int n=1;n<=8;n++){
      int nx=x+rd::d8x[n], ny=y+rd::d8y[n];
      if(!dem.inGrid(nx,ny)) continue;
      if(sd(x,y)!=NO_FLOW && sd(nx,ny)==NO_FLOW && dem(nx,ny)==dem(x,y)){ low_edges.push_back({x,y}); break; }
      if(sd(x,y)==NO_FLOW && dem(x,y)<dem(nx,ny)){ high_edges.push_back({x,y}); break; }
    }

  const auto same_flat = [&](int x,int y,int nx,int ny){          // == "same flat label"
    return dem.inGrid(nx,ny) && sd(nx,ny)==NO_FLOW && dem(nx,ny)==dem(x,y);
  };

  // away_dist: multi-source min-BFS from high edges (richdem's BuildAwayGradient, label->elevation).
  rd::Array2D<int32_t> away(W,H,0);
  { std::deque<rd::GridCell> q=high_edges; q.push_back({-1,-1}); int loops=1;
    while(q.size()!=1){ auto c=q.front(); q.pop_front();
      if(c.x==-1){ loops++; q.push_back({-1,-1}); continue; }
      if(away(c.x,c.y)>0) continue; away(c.x,c.y)=loops;
      for(int n=1;n<=8;n++){ int nx=c.x+rd::d8x[n], ny=c.y+rd::d8y[n]; if(same_flat(c.x,c.y,nx,ny)) q.push_back({nx,ny}); } } }

  // flat_height: max-relaxation of away over same-elevation adjacency (label-free).
  rd::Array2D<int32_t> fh(W,H,0);
  for(int i=0;i<(int)fh.size();i++) fh(i)=away(i);
  for(bool changed=true; changed; ){ changed=false;
    for(int x=0;x<W;x++) for(int y=0;y<H;y++){ if(dem.isNoData(x,y)) continue;
      for(int n=1;n<=8;n++){ int nx=x+rd::d8x[n], ny=y+rd::d8y[n];
        if(dem.inGrid(nx,ny) && dem(nx,ny)==dem(x,y) && fh(nx,ny)>fh(x,y)){ fh(x,y)=fh(nx,ny); changed=true; } } } }

  // towards + combine: mask = flat_height - away + 2*loops (richdem's BuildTowardsCombinedGradient).
  rd::Array2D<int32_t> mask(W,H,0);
  for(int i=0;i<(int)mask.size();i++) mask(i) = -away(i);
  { std::deque<rd::GridCell> q=low_edges; q.push_back({-1,-1}); int loops=1;
    while(q.size()!=1){ auto c=q.front(); q.pop_front();
      if(c.x==-1){ loops++; q.push_back({-1,-1}); continue; }
      if(mask(c.x,c.y)>0) continue;
      mask(c.x,c.y) = (mask(c.x,c.y)!=0) ? (fh(c.x,c.y)+mask(c.x,c.y)+2*loops) : (2*loops);
      for(int n=1;n<=8;n++){ int nx=c.x+rd::d8x[n], ny=c.y+rd::d8y[n]; if(same_flat(c.x,c.y,nx,ny)) q.push_back({nx,ny}); } } }

  // (A) == (B) on flat cells; also confirm the label-free flat_height matches the label-based one.
  long flats=0, mask_diff=0, fh_diff=0;
  std::vector<int32_t> fh_label(lab.max()+1, 0);
  for(int i=0;i<(int)away.size();i++) if(sd(i)==NO_FLOW && lab(i)>0) fh_label[lab(i)] = std::max(fh_label[lab(i)], away(i));
  for(int i=0;i<(int)fm_ref.size();i++){
    if(sd(i)!=NO_FLOW) continue; flats++;
    if(mask(i)!=fm_ref(i)) mask_diff++;
    if(lab(i)>0 && fh(i)!=fh_label[lab(i)]) fh_diff++;
  }
  // Also check the end-to-end flowdir producer: resolve_flat_flowdirs_option2 (the label-free
  // relaxation form used by ENH-1) must equal the reference resolve_flat_flowdirs, flowdir for flowdir.
  rd::Array2D<int8_t> fa=sd, fb=sd;
  resolve_flat_flowdirs(dem, fa);            // richdem-based reference (full grid)
  resolve_flat_flowdirs_option2(dem, fb);    // ENH-1 relaxation reconstruction
  long fd_diff=0; for(int i=0;i<(int)fa.size();i++) if(fa(i)!=fb(i)) fd_diff++;

  const bool ok = (mask_diff==0 && fh_diff==0 && fd_diff==0);
  std::cout<<(ok?"FLAT-RECONSTRUCT-OK":"FLAT-RECONSTRUCT-DIFFER")
           <<" flat_cells="<<flats<<" high_edges="<<high_edges.size()<<" low_edges="<<low_edges.size()
           <<" mask_diff="<<mask_diff<<" label_free_flat_height_diff="<<fh_diff
           <<" option2_flowdir_diff="<<fd_diff<<"\n";
  return ok?0:1;
}
