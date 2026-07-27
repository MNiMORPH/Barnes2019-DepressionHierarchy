// dephier_mpi -- distributed DepressionHierarchy on the thread-comm shim (comm_thread.hpp),
// the no-cluster validation transport that ports to real MPI unchanged. Built incrementally,
// each increment diffed against an in-process oracle exactly as tools/dephier_stitch.cpp was.
//
// Design principle (Wickert): keep everything LOCAL as much as possible -- each rank owns one
// tile, runs on its own columns, and only 1-cell perimeter strips cross a seam. Nothing gathers
// a per-cell field to rank 0; only the O(boundary) linkage and (later) the meta-tree do. This
// mirrors richdem's parallel_priority_flood, which writes each filled tile to disk itself and
// ships only its edge/graph up (main.cpp:328 saveGDAL vs :654 CommSend(job1)).
//
// INCREMENT 1 (this file): per-rank Phase A/B with a halo-exchanged BOUNDARY pre-label.
//   * component 2 (per-rank PhaseAB) + component 3 (perimeter-strip exchange) of the plan.
//   * Each rank extracts its own DEM columns, exchanges its 1-column edge elevations with its
//     seam neighbours via the shim, computes the BOUNDARY pre-label from own+halo columns
//     (no full grid), and runs GetDepressionHierarchyPhaseAB locally.
//   * ORACLE: the same per-tile PhaseAB but with BOUNDARY computed from the FULL grid (the
//     tools/dephier_stitch.cpp per-tile step). The distributed label+fd must be bit-identical.
//   Verifies: (a) the shim moves DH edge strips across a multi-rank line with no deadlock, and
//   (b) BOUNDARY-from-a-1-column-halo == BOUNDARY-from-the-full-grid (the local-first claim).

#include "tools/comm_thread.hpp"

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;
namespace c  = commt;

using dh::dh_label_t;
using dh::OCEAN;
using dh::NO_DEP;
using dh::BOUNDARY;

// ---- small helpers (shared with the in-process stitch; kept identical) ----
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

// A 1-column perimeter strip: elevation + ocean flag down a seam-edge column. This is the
// message that crosses a rank boundary (cereal-serialized by the shim). Ocean is shipped
// explicitly so the receiver's is_ocean() is bit-identical to the sender's (no re-derivation
// from a possibly-NaN noData sentinel).
struct EdgeStrip {
  std::vector<float>   elev;
  std::vector<uint8_t> ocean;
  template<class Ar> void serialize(Ar &ar){ ar(elev, ocean); }
};

int main(int argc, char **argv){
  if(argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level> <Split Cols (comma-sep)>\n";
    return -1;
  }
  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);

  rd::Array2D<float> full(in_name);
  const int W = full.width(), H = full.height();

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
  bounds.push_back(W);
  const int ntiles = bounds.size() - 1;

  // Split columns must be strictly increasing and interior to the grid, else a tile would
  // extend past the DEM (out-of-grid read). Fail cleanly rather than abort.
  for(size_t i=1;i+1<bounds.size();i++)
    if(bounds[i]<=bounds[i-1] || bounds[i]>=W){
      std::cerr<<"error: split column "<<bounds[i]<<" is out of range (0 < col < width="<<W
               <<") or not strictly increasing\n";
      return -1;
    }

  const auto tile_of = [&](int x){ int t=0; while(t+1<(int)bounds.size() && x>=bounds[t+1]) t++; return t; };
  const auto is_ocean_full = [&](int x,int y){ return full.isNoData(x,y) || full(x,y)==ocean_level; };

  // ---- ORACLE: per-tile PhaseAB with BOUNDARY from the FULL grid (dephier_stitch's step) ----
  // Lowest strictly-downhill land neighbour on the full grid; ties by highest cell index
  // (radix pop order). Identical to the stitch's drain().
  const auto drain_full = [&](int x,int y,int &lx,int &ly,bool &to_ocean)->bool{
    const float focal = full(x,y);
    float best = std::numeric_limits<float>::infinity(); bool found=false; to_ocean=false;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
      if(!dx && !dy) continue;
      int nx=x+dx, ny=y+dy;
      if(!full.inGrid(nx,ny)) continue;
      if(is_ocean_full(nx,ny)){ to_ocean=true; continue; }
      const float e = full(nx,ny);
      if(e >= focal) continue;
      if(e <= best){ best=e; lx=nx; ly=ny; found=true; }
    }
    return found;
  };

  struct Result {
    rd::Array2D<dh_label_t> label;      // tile-LOCAL labels (OCEAN/BOUNDARY sentinels + local deps)
    rd::Array2D<int8_t>     fd;
    rd::Array2D<dh_label_t> glab;        // tile slice remapped into the GLOBAL namespace (pre-conduit)
    int        nboundary = 0;
    int        ndep      = 0;            // local depressions, excl. ocean node 0 (deps.size()-1)
    dh_label_t offset    = 0;            // this tile's global label offset (prefix sum)
  };
  std::vector<Result> oracle(ntiles), dist(ntiles);

  // Local label -> global label. Sentinels (OCEAN, BOUNDARY) are global as-is; a real local
  // depression k maps to offset+(k-1). Identical to Tile::g() in the stitch.
  const auto gmap = [&](dh_label_t local, dh_label_t offset)->dh_label_t{
    if(local==OCEAN || local==BOUNDARY) return local;
    return offset + (local - 1);
  };

  for(int t=0;t<ntiles;t++){
    const int x0=bounds[t], x1=bounds[t+1];
    rd::Array2D<float>     dem   = extract_cols(full, x0, x1);
    rd::Array2D<dh_label_t> label = ocean_labels(dem, ocean_level);
    std::vector<int> seam_cols;
    if(t>0)        seam_cols.push_back(x0);
    if(t<ntiles-1) seam_cols.push_back(x1-1);
    int nb=0;
    for(int gc : seam_cols)
      for(int y=0;y<H;y++){
        const int lc = gc - x0;
        if(label(lc,y)!=NO_DEP) continue;
        int lx,ly; bool to_ocean;
        if(drain_full(gc,y,lx,ly,to_ocean) && !to_ocean && tile_of(lx)!=t){ label(lc,y)=BOUNDARY; nb++; }
      }
    rd::Array2D<int8_t> fd(dem.width(), dem.height(), rd::NO_FLOW);
    dh::DepressionHierarchy<float> deps; std::vector<dh::Outlet<float>> outlets;
    dh::GetDepressionHierarchyPhaseAB<float,rd::Topology::D8>(dem, label, fd, deps, outlets);
    oracle[t].label = std::move(label);
    oracle[t].fd    = std::move(fd);
    oracle[t].nboundary = nb;
    oracle[t].ndep  = (int)deps.size() - 1;
  }
  // Oracle namespace remap: prefix-sum the per-tile depression counts (the stitch's next_offset
  // accumulation), then map each tile's local labels into the global grid (BOUNDARY unresolved).
  rd::Array2D<dh_label_t> gLabel_oracle(W, H, OCEAN);
  {
    dh_label_t next = 1;
    for(int t=0;t<ntiles;t++){ oracle[t].offset = next; next += oracle[t].ndep; }
    for(int t=0;t<ntiles;t++){
      const int x0=bounds[t], x1=bounds[t+1];
      for(int y=0;y<H;y++) for(int x=x0;x<x1;x++)
        gLabel_oracle(x,y) = gmap(oracle[t].label(x-x0,y), oracle[t].offset);
    }
  }

  // ---- DISTRIBUTED: one thread-rank per tile, BOUNDARY from a 1-column halo ----
  enum { TAG_L2R=1, TAG_R2L=2, TAG_COUNT=3, TAG_OFFSET=4 };   // strip dirs; count->rank0; offset->rank r
  auto rank_main = [&](){
    const int r  = c::CommRank();
    const int x0 = bounds[r], x1 = bounds[r+1], w = x1-x0;
    rd::Array2D<float>      dem   = extract_cols(full, x0, x1);   // this rank's own columns only
    rd::Array2D<dh_label_t> label = ocean_labels(dem, ocean_level);

    // Build this rank's edge strips and exchange with seam neighbours (non-blocking send
    // into the neighbour's inbox, then blocking recv -- no deadlock).
    const auto strip_of = [&](int lc){
      EdgeStrip s; s.elev.resize(H); s.ocean.resize(H);
      for(int y=0;y<H;y++){ s.elev[y]=dem(lc,y); s.ocean[y]= (dem.isNoData(lc,y)||dem(lc,y)==ocean_level); }
      return s;
    };
    if(r+1<ntiles) c::CommSend(strip_of(w-1), r+1, TAG_L2R);      // my right edge -> r+1's left halo
    if(r-1>=0)     c::CommSend(strip_of(0),   r-1, TAG_R2L);      // my left  edge -> r-1's right halo
    EdgeStrip haloL, haloR;
    if(r-1>=0)     c::CommRecv(haloL, r-1, TAG_L2R);              // column x0-1
    if(r+1<ntiles) c::CommRecv(haloR, r+1, TAG_R2L);              // column x1

    // Local elevation/ocean access over [x0-1, x1] (own columns + the two halos).
    const auto elev = [&](int gx,int y)->float{
      if(gx>=x0 && gx<x1) return dem(gx-x0,y);
      if(gx==x0-1)        return haloL.elev[y];
      return haloR.elev[y];                                       // gx==x1
    };
    const auto is_ocean = [&](int gx,int y)->bool{
      if(gx>=x0 && gx<x1) return dem.isNoData(gx-x0,y) || dem(gx-x0,y)==ocean_level;
      if(gx==x0-1)        return haloL.ocean[y];
      return haloR.ocean[y];
    };
    // Rank-local drain: same rule as drain_full, but reading only own columns + halos.
    const auto drain_local = [&](int gx,int y,int &lx,bool &to_ocean)->bool{
      const float focal = elev(gx,y);
      float best = std::numeric_limits<float>::infinity(); bool found=false; to_ocean=false;
      for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
        if(!dx && !dy) continue;
        const int nx=gx+dx, ny=y+dy;
        if(nx<0||nx>=W||ny<0||ny>=H) continue;
        if(is_ocean(nx,ny)){ to_ocean=true; continue; }
        const float e = elev(nx,ny);
        if(e >= focal) continue;
        if(e <= best){ best=e; lx=nx; found=true; }              // <= : higher-index tie wins
      }
      return found;
    };

    std::vector<int> seam_cols;
    if(r>0)        seam_cols.push_back(x0);
    if(r<ntiles-1) seam_cols.push_back(x1-1);
    int nb=0;
    for(int gc : seam_cols)
      for(int y=0;y<H;y++){
        const int lc = gc - x0;
        if(label(lc,y)!=NO_DEP) continue;
        int lx; bool to_ocean;
        if(drain_local(gc,y,lx,to_ocean) && !to_ocean && tile_of(lx)!=r){ label(lc,y)=BOUNDARY; nb++; }
      }

    rd::Array2D<int8_t> fd(dem.width(), dem.height(), rd::NO_FLOW);
    dh::DepressionHierarchy<float> deps; std::vector<dh::Outlet<float>> outlets;
    dh::GetDepressionHierarchyPhaseAB<float,rd::Topology::D8>(dem, label, fd, deps, outlets);

    // Namespace remap (eng-doc component 6): gather per-tile depression counts to rank 0,
    // prefix-sum into global offsets, scatter each rank its offset. This is the shim analogue
    // of MPI_Allgather(count) + prefix sum -- O(ntiles) tiny messages, no per-cell data.
    const int mycount = (int)deps.size() - 1;
    dh_label_t myoffset = 1;
    if(r==0){
      std::vector<int> counts(ntiles); counts[0]=mycount;
      for(int t=1;t<ntiles;t++) c::CommRecv(counts[t], t, TAG_COUNT);
      dh_label_t off = 1;
      for(int t=0;t<ntiles;t++){                                  // prefix sum in tile order
        const dh_label_t this_off = off; off += counts[t];
        if(t==0) myoffset = this_off; else c::CommSend(this_off, t, TAG_OFFSET);
      }
    } else {
      c::CommSend(mycount, 0, TAG_COUNT);
      c::CommRecv(myoffset, 0, TAG_OFFSET);
    }

    // Map this rank's local labels into the global namespace (its own slice; no gather).
    rd::Array2D<dh_label_t> glab(w, H);
    for(int y=0;y<H;y++) for(int x=0;x<w;x++) glab(x,y) = gmap(label(x,y), myoffset);

    dist[r].label = std::move(label);                            // ranks write disjoint slots
    dist[r].fd    = std::move(fd);
    dist[r].glab  = std::move(glab);
    dist[r].nboundary = nb;
    dist[r].ndep  = mycount;
    dist[r].offset= myoffset;
    c::CommBarrier();
  };
  c::CommInit(ntiles, rank_main);

  // ---- verify: distributed label+fd must equal the full-grid oracle, cell for cell ----
  long ldiff=0, fdiff=0, cells=0, nb_o=0, nb_d=0;
  for(int t=0;t<ntiles;t++){
    nb_o += oracle[t].nboundary; nb_d += dist[t].nboundary;
    for(unsigned i=0;i<oracle[t].label.size();i++){
      cells++;
      if(oracle[t].label(i)!=dist[t].label(i)) ldiff++;
      if(oracle[t].fd(i)   !=dist[t].fd(i))    fdiff++;
    }
  }
  const bool phaseab_ok = (ldiff==0 && fdiff==0);
  std::cout<<(phaseab_ok ? "MPI-PHASEAB-MATCH " : "MPI-PHASEAB-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" label_diff="<<ldiff<<" fd_diff="<<fdiff<<"/"<<cells
           <<" boundary(oracle="<<nb_o<<" dist="<<nb_d<<")\n";

  // ---- remap check: global offsets + the assembled global label grid (pre-conduit) must
  // equal the stitch's assembly, cell for cell. ----
  long odiff=0, gdiff=0;
  rd::Array2D<dh_label_t> gLabel_dist(W, H, OCEAN);
  for(int t=0;t<ntiles;t++){
    if(dist[t].offset != oracle[t].offset) odiff++;
    const int x0=bounds[t], x1=bounds[t+1];
    for(int y=0;y<H;y++) for(int x=x0;x<x1;x++)
      gLabel_dist(x,y) = dist[t].glab(x-x0,y);
  }
  for(unsigned i=0;i<gLabel_oracle.size();i++)
    if(gLabel_oracle(i)!=gLabel_dist(i)) gdiff++;
  const bool remap_ok = (odiff==0 && gdiff==0);
  std::cout<<(remap_ok ? "MPI-REMAP-MATCH " : "MPI-REMAP-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" offset_diff="<<odiff<<"/"<<ntiles
           <<" glabel_diff="<<gdiff<<"/"<<cells<<" n_global="<<(oracle[ntiles-1].offset+oracle[ntiles-1].ndep)<<"\n";

  return (phaseab_ok && remap_ok) ? 0 : 1;
}
