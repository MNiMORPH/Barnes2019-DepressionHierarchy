// dh_flats -- Barnes-2014 flat-flowdir resolution, shared by the in-process stitch
// (tools/dephier_stitch.cpp) and the distributed harness (tools/dephier_mpi.cpp). Extracted
// verbatim from the stitch; no behaviour change. See PARALLEL_DEPHIER_ENGINEERING.md section 6.
//
// KEY EQUIVALENCE for the distributed flat resolution (ENH-1, ENHANCEMENTS.md): richdem's
// resolve_flats_barnes labels flats with label_this, a flood-fill over equal-elevation
// D8-neighbours -- so "same flat label" is IDENTICAL to "D8-adjacent AND same elevation." That
// lets the whole flat_mask be rebuilt with NO global flat labeling, as three order-independent
// relaxations over same-elevation adjacency (away-BFS from high edges, towards-BFS from low
// edges, and a max-relaxation for flat_height), each of which distributes via a 1-column seam
// exchange iterated to convergence. Proven bit-identical to resolve_flats_barnes by the
// executable regression test tools/flat_mask_reconstruct_test.cpp (CTest flat_mask_reconstruct).
#pragma once

#include <richdem/common/Array2D.hpp>
#include <richdem/flowmet/d8_flowdirs.hpp>
#include <richdem/flats/flat_resolution.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace rd = richdem;

// Overlay Barnes-2014 resolve_flats directions onto the FLAT cells of `fd` (cells with
// no strictly-lower neighbour), leaving non-flat cells untouched. DH assigns flat cells
// an order-dependent flood-claim direction; this replaces it with the geometry-
// deterministic flat routing (BuildAwayGradient + BuildTowardsCombinedGradient -> mask
// -> d8_flow_flats), which is reproducible and therefore agrees serial vs. distributed.
// See PARALLEL_DEPHIER_ENGINEERING.md section 6. (This is the full-grid form; the
// distributed form computes the two gradient BFSs with a cross-seam boundary exchange.)
inline void resolve_flat_flowdirs(const rd::Array2D<float> &dem, rd::Array2D<int8_t> &fd){
  rd::Array2D<int8_t> sd(dem.width(), dem.height(), rd::NO_FLOW);
  rd::d8_flow_directions(dem, sd);                    // steepest descent: NO_FLOW on flats+pits
  std::vector<char> was_flat(dem.size(), 0);
  for(unsigned i=0;i<dem.size();i++)
    if(!dem.isNoData(i) && sd(i)==rd::NO_FLOW) was_flat[i]=1;   // flat or genuine pit
  rd::Array2D<int32_t> flat_mask, labels;
  rd::resolve_flats_barnes(dem, sd, flat_mask, labels);
  rd::d8_flow_flats(flat_mask, labels, sd);           // fill flat directions from the mask
  for(unsigned i=0;i<dem.size();i++)
    if(was_flat[i]) fd(i)=sd(i);                      // pits stay NO_FLOW; flats get resolved dir
}

// MOVE 2: footprint-bounded distributed flat resolution. Instead of resolving flats on
// the full grid, each tile resolves its OWN flats with resolve_flats on the tile plus an
// adaptive boundary halo, grown until the owned region stops changing -- a purely local
// convergence test (no full grid). This reproduces the full-grid resolve_flats result
// exactly (validated bit-identical), because the geometry-deterministic routing is
// locally determined away from the seam; the halo need only reach the flat's cross-seam
// extent. Overlays the resolved directions onto the flat cells of `fd`, like the
// full-grid form. See PARALLEL_DEPHIER_ENGINEERING.md section 6.
inline constexpr int8_t NOT_FLAT = -2;   // sentinel: cell is not a flat (leave its flowdir alone)
inline void resolve_flat_flowdirs_into(const rd::Array2D<float> &sub, rd::Array2D<int8_t> &out){
  rd::Array2D<int8_t> sd(sub.width(), sub.height(), rd::NO_FLOW);
  rd::d8_flow_directions(sub, sd);
  std::vector<char> flat(sub.size(), 0);
  for(unsigned i=0;i<sub.size();i++) if(!sub.isNoData(i) && sd(i)==rd::NO_FLOW) flat[i]=1;
  rd::Array2D<int32_t> fm,lb; rd::resolve_flats_barnes(sub,sd,fm,lb); rd::d8_flow_flats(fm,lb,sd);
  out = rd::Array2D<int8_t>(sub.width(), sub.height(), NOT_FLAT);
  for(unsigned i=0;i<sub.size();i++) if(flat[i]) out(i)=sd(i);   // flat cell -> resolved dir (0=mesa..8)
}
// `halo_cap` bounds the halo (option 3, PARALLEL_DEPHIER_ENGINEERING.md section 6): the
// halo grows until the owned region stabilizes (bit-identical) OR the cap, whichever comes
// first. A flat WIDER than the cap is resolved with a capped halo -> a valid convergent
// flow field (same tree, same sinks) but not necessarily serial-identical in that flat's
// interior. This guarantees per-rank footprint O(N/P)+O(cap*boundary) for ANY input, so
// nothing can break the build; the bit-identical option-2 exchange is a later optimization.
// Returns the number of tiles whose halo hit the cap without converging. A cap >= grid
// width is effectively unlimited (grows to the full grid = exact); that is the default.
inline int resolve_flat_flowdirs_distributed(const rd::Array2D<float> &full,
                                             const std::vector<int> &bounds,
                                             rd::Array2D<int8_t> &fd,
                                             int halo_cap){
  const int W=full.width(), H=full.height();
  const auto owned_halo = [&](int x0,int x1,int h,rd::Array2D<int8_t> &owned){ // owned flat dirs, halo h
    const int a=std::max(0,x0-h), b=std::min(W,x1+h), w=b-a;
    rd::Array2D<float> sub(w,H); sub.setNoData(full.noData());
    for(int y=0;y<H;y++) for(int x=a;x<b;x++) sub(x-a,y)=full(x,y);
    rd::Array2D<int8_t> sfd; resolve_flat_flowdirs_into(sub, sfd);
    owned = rd::Array2D<int8_t>(W,H,NOT_FLAT);
    for(int y=0;y<H;y++) for(int x=x0;x<x1;x++) owned(x,y)=sfd(x-a,y);
  };
  int capped=0;
  for(size_t t=0;t+1<bounds.size();t++){
    const int x0=bounds[t], x1=bounds[t+1];
    const int hmax=std::min(halo_cap, W);
    rd::Array2D<int8_t> prev, cur;
    int h=std::min(2,hmax); owned_halo(x0,x1,h,prev);
    bool converged=false;
    while(h<hmax){                                      // grow until owned region stable OR the cap
      const int h2=std::min(h*2,hmax);
      owned_halo(x0,x1,h2,cur);
      bool same=true;
      for(int y=0;y<H && same;y++) for(int x=x0;x<x1;x++) if(prev(x,y)!=cur(x,y)){ same=false; break; }
      prev=cur; h=h2;
      if(same){ converged=true; break; }
    }
    if(!converged) capped++;                            // a flat wider than the cap: valid, maybe not identical
    for(int y=0;y<H;y++) for(int x=x0;x<x1;x++)          // overlay resolved flat directions
      if(prev(x,y)!=NOT_FLAT) fd(x,y)=prev(x,y);         // incl. NO_FLOW (mesa), matching full-grid
  }
  return capped;
}
