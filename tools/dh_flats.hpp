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

// ---- ENH-1 / "option 2": label-free flat resolution as three order-independent RELAXATIONS over
// same-elevation adjacency (see the equivalence note at the top of this file; proven bit-identical by
// tools/flat_mask_reconstruct_test.cpp). Setup (`flat_seed`) and finish (`flat_finish`: mask ->
// masked flowdir) are shared; only the relaxation SCHEDULING differs between the whole-grid form and
// the per-tile + 1-column-seam-exchange form. Both overlay onto the flat cells of `fd` exactly like
// resolve_flat_flowdirs.
inline bool flat_at(const rd::Array2D<float>&dem, const rd::Array2D<int8_t>&sd, int x,int y){
  return sd.inGrid(x,y) && !dem.isNoData(x,y) && sd(x,y)==rd::NO_FLOW;   // flat/pit cell (NO_FLOW)
}
// Fill sd (d8) and the relaxation seeds. HIGH edge (flat cell touching higher terrain) -> away=1;
// LOW edge (non-flat cell touching a same-elevation flat) -> tw=1; everything else INF.
inline void flat_seed(const rd::Array2D<float>&dem, rd::Array2D<int8_t>&sd,
                      rd::Array2D<int32_t>&away, rd::Array2D<int32_t>&tw, int32_t INF){
  const int W=dem.width(), H=dem.height();
  rd::d8_flow_directions(dem, sd);
  away.setAll(INF); tw.setAll(INF);
  for(int y=0;y<H;y++) for(int x=0;x<W;x++){
    if(dem.isNoData(x,y)) continue;
    for(int n=1;n<=8;n++){ int nx=x+rd::d8x[n], ny=y+rd::d8y[n]; if(!dem.inGrid(nx,ny)||dem.isNoData(nx,ny)) continue;
      if(flat_at(dem,sd,x,y)  && dem(nx,ny)>dem(x,y))                             away(x,y)=1;
      if(!flat_at(dem,sd,x,y) && flat_at(dem,sd,nx,ny) && dem(nx,ny)==dem(x,y))   tw(x,y)=1;
    }
  }
}
// Mask (Barnes-2014) then masked steepest-descent, overlaid onto EVERY flat cell (fd(was_flat)=dir,
// NO_FLOW where border/unresolved), so the result is independent of what `fd` held on flats before.
inline void flat_finish(const rd::Array2D<float>&dem, const rd::Array2D<int8_t>&sd,
                        const rd::Array2D<int32_t>&away, const rd::Array2D<int32_t>&tw,
                        const rd::Array2D<int32_t>&fh, rd::Array2D<int8_t>&fd, int32_t INF){
  const int W=dem.width(), H=dem.height();
  rd::Array2D<int32_t> mask(W,H,INF);
  for(int y=0;y<H;y++) for(int x=0;x<W;x++){
    const bool f=flat_at(dem,sd,x,y);
    if(f && tw(x,y)<INF)       mask(x,y) = (away(x,y)<INF)? (fh(x,y)-away(x,y)+2*tw(x,y)) : (2*tw(x,y));
    else if(!f && tw(x,y)<INF) mask(x,y) = 2*tw(x,y);                             // low edge (=2)
  }
  rd::Array2D<int8_t> out(W,H,rd::NO_FLOW);
  for(int y=1;y<H-1;y++) for(int x=1;x<W-1;x++){                                  // interior, like d8_flow_flats
    if(!flat_at(dem,sd,x,y) || mask(x,y)>=INF) continue;
    int best=mask(x,y), dir=rd::NO_FLOW;
    for(int n=1;n<=8;n++){ int nx=x+rd::d8x[n], ny=y+rd::d8y[n];
      if(dem(nx,ny)!=dem(x,y) || mask(nx,ny)>=INF) continue;                      // same-elevation (== same label)
      if(mask(nx,ny)<best || (mask(nx,ny)==best && dir>0 && dir%2==0 && n%2==1)){ best=mask(nx,ny); dir=n; } }
    out(x,y)=(int8_t)dir;
  }
  for(int y=0;y<H;y++) for(int x=0;x<W;x++) if(flat_at(dem,sd,x,y)) fd(x,y)=out(x,y);
}

// WHOLE-GRID form: each relaxation sweeps the whole grid to its fixed point (the distributed algorithm
// scheduled globally). Used to validate bit-identity and as the in-process harness flat pass.
inline void resolve_flat_flowdirs_option2(const rd::Array2D<float>&dem, rd::Array2D<int8_t>&fd){
  const int W=dem.width(), H=dem.height(); const int32_t INF=1<<28;
  rd::Array2D<int8_t> sd(W,H,rd::NO_FLOW);
  rd::Array2D<int32_t> away(W,H,INF), tw(W,H,INF), fh(W,H,0);
  flat_seed(dem,sd,away,tw,INF);
  const auto rmin=[&](rd::Array2D<int32_t>&d){ for(bool go=true; go; ){ go=false;
    for(int y=0;y<H;y++) for(int x=0;x<W;x++){ if(!flat_at(dem,sd,x,y)) continue;
      for(int n=1;n<=8;n++){ int nx=x+rd::d8x[n], ny=y+rd::d8y[n];
        if(dem.inGrid(nx,ny)&&!dem.isNoData(nx,ny)&&dem(nx,ny)==dem(x,y)&&d(nx,ny)+1<d(x,y)){ d(x,y)=d(nx,ny)+1; go=true; } } } } };
  rmin(away); rmin(tw);
  for(int y=0;y<H;y++) for(int x=0;x<W;x++) fh(x,y)=(flat_at(dem,sd,x,y)&&away(x,y)<INF)?away(x,y):0;
  for(bool go=true; go; ){ go=false;
    for(int y=0;y<H;y++) for(int x=0;x<W;x++){ if(dem.isNoData(x,y)) continue;
      for(int n=1;n<=8;n++){ int nx=x+rd::d8x[n], ny=y+rd::d8y[n];
        if(dem.inGrid(nx,ny)&&!dem.isNoData(nx,ny)&&dem(nx,ny)==dem(x,y)&&fh(nx,ny)>fh(x,y)){ fh(x,y)=fh(nx,ny); go=true; } } } }
  flat_finish(dem,sd,away,tw,fh,fd,INF);
}

// PER-TILE + 1-column-seam-exchange form (ENH-1's genuine distributed realization). Each ROUND every
// tile relaxes its OWNED columns [x0,x1) to local convergence, reading cross-seam neighbours from the
// PREVIOUS round's snapshot -- so a value crosses at most one seam per round; rounds ~ flat diameter
// in tiles, and only O(boundary) (the seam columns) crosses per round. No tile reads another tile's
// interior. Bit-identical to the whole-grid form (same monotone min/max fixed point). Returns the
// total relaxation rounds (a diagnostic: the seam-exchange count). `bounds` are the tile column splits.
inline int resolve_flat_flowdirs_option2_tiled(const rd::Array2D<float>&dem, const std::vector<int>&bounds,
                                               rd::Array2D<int8_t>&fd){
  const int W=dem.width(), H=dem.height(); const int32_t INF=1<<28; const int nt=(int)bounds.size()-1;
  rd::Array2D<int8_t> sd(W,H,rd::NO_FLOW);
  rd::Array2D<int32_t> away(W,H,INF), tw(W,H,INF), fh(W,H,0);
  flat_seed(dem,sd,away,tw,INF);
  int rounds=0;
  const auto tiled=[&](rd::Array2D<int32_t>&d, bool is_max, bool flat_only){
    for(bool go=true; go; ){ go=false; rounds++;
      rd::Array2D<int32_t> prev=d;                             // snapshot: cross-seam reads see last round
      for(int t=0;t<nt;t++){ const int x0=bounds[t], x1=bounds[t+1];
        for(bool lgo=true; lgo; ){ lgo=false;                  // local convergence within this tile
          for(int y=0;y<H;y++) for(int x=x0;x<x1;x++){
            if(flat_only ? !flat_at(dem,sd,x,y) : dem.isNoData(x,y)) continue;
            for(int n=1;n<=8;n++){ int nx=x+rd::d8x[n], ny=y+rd::d8y[n];
              if(!dem.inGrid(nx,ny)||dem.isNoData(nx,ny)||dem(nx,ny)!=dem(x,y)) continue;
              const int32_t nv = (nx>=x0 && nx<x1) ? d(nx,ny) : prev(nx,ny);   // in-tile live / cross-seam snapshot
              if(is_max){ if(nv>d(x,y)){ d(x,y)=nv;   lgo=true; } }
              else      { if(nv+1<d(x,y)){ d(x,y)=nv+1; lgo=true; } }
            }
          }
        }
      }
      for(int i=0;i<(int)d.size();i++) if(d(i)!=prev(i)){ go=true; break; }
    }
  };
  tiled(away,false,true);
  tiled(tw,  false,true);
  for(int y=0;y<H;y++) for(int x=0;x<W;x++) fh(x,y)=(flat_at(dem,sd,x,y)&&away(x,y)<INF)?away(x,y):0;
  tiled(fh, true, false);
  flat_finish(dem,sd,away,tw,fh,fd,INF);
  return rounds;
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
