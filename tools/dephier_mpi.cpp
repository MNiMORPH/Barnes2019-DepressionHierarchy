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

// Comm backend: the thread shim (in-process validation) by default, or real MPI when built with
// -DDH_USE_MPI (one rank per process; launch with `mpirun -np <ntiles>`). Both expose the same
// commt API, so the algorithm code below is identical -- only the transport and driver differ.
#ifdef DH_USE_MPI
#include "tools/comm_mpi.hpp"
#else
#include "tools/comm_thread.hpp"
#endif
#include "dh_canonical.hpp"   // dhtest::canonicalize / invariants (shared with the stitch)
#include "dh_collapse.hpp"    // CollapseSeamArtifacts            (shared with the stitch)
#include "dh_flats.hpp"       // resolve_flat_flowdirs[_*]        (shared with the stitch)

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;
namespace c  = commt;

// Non-intrusive cereal serializer for a Depression, so leaf records can cross the shim to
// rank 0 (eng-doc component 7). Found by ADL (Depression lives in richdem::dephier); does NOT
// modify the published library. Carries every field so the gathered record is a faithful copy.
namespace richdem { namespace dephier {
template<class Ar, class E>
void serialize(Ar &ar, Depression<E> &d){
  ar(d.pit_cell, d.out_cell, d.parent, d.odep, d.geolink, d.pit_elev, d.out_elev,
     d.lchild, d.rchild, d.ocean_parent, d.ocean_linked, d.dep_label,
     d.cell_count, d.dep_vol, d.water_vol, d.total_elevation);
}
}}

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
  std::vector<uint8_t> nodata;   // distinguish a noData cell from an ocean-valued land neighbour
  template<class Ar> void serialize(Ar &ar){ ar(elev, ocean, nodata); }
};

// Conduit resolution messages (eng-doc section 2), all keyed on a global (gx,gy) cell.
//   LWRec  -- a BOUNDARY cell's Phase-1 result: a terminal GLOBAL label, or an EXIT to a
//             cross-seam cell (ex,ey). One per BOUNDARY cell (O(boundary) per rank).
//   CellLab-- (gx,gy)->global label; used both for shipping a rank's edge-column labels
//             (so rank 0 can read an EXIT target's label) and for scattering resolved labels.
struct LWRec  { int32_t gx,gy; uint8_t terminal; dh_label_t label; int32_t ex,ey;
                template<class Ar> void serialize(Ar &ar){ ar(gx,gy,terminal,label,ex,ey); } };
struct CellLab{ int32_t gx,gy; dh_label_t label;
                template<class Ar> void serialize(Ar &ar){ ar(gx,gy,label); } };

// Outlet-set messages (eng-doc component 4).
//   ResStrip -- a resolved edge column {global label, elevation, noData flag}, exchanged so a
//               rank can run Barnes' HandleEdge across its seam without the neighbour's interior.
//   ORec     -- one reduced outlet-DB entry: a depression pair + the outlet (higher-of-pair
//               elevation and its cell). O(#local depression-adjacencies) per rank, not O(N/P).
struct ResStrip{ std::vector<dh_label_t> label; std::vector<float> elev; std::vector<uint8_t> nodata;
                 template<class Ar> void serialize(Ar &ar){ ar(label,elev,nodata); } };
struct ORec    { dh_label_t depa, depb; float oelev; int64_t ocell;
                 template<class Ar> void serialize(Ar &ar){ ar(depa,depb,oelev,ocell); } };

// A rank's per-cell grid slices, gathered to rank 0 under real MPI so the FULL per-cell verification
// (not just the tree) can run there. Array2D has no cereal serializer, so ship the raw row-major data
// + dims. O(N) to rank 0 -- validation-grade only (the algorithm itself never gathers a per-cell field).
struct DistSlice {
  int32_t w=0, h=0; int32_t nboundary=0, ndep=0; dh_label_t offset=0;
  std::vector<dh_label_t> label, glab, glab_pc;
  std::vector<int8_t>     fd, gfix;
  template<class Ar> void serialize(Ar &ar){ ar(w,h,nboundary,ndep,offset,label,glab,glab_pc,fd,gfix); }
};

// Seam exchange + convergence for the per-rank flat resolution (dh_flats.hpp resolve_flat_flowdirs_rank).
// exch(): send my owned edge columns to my seam neighbours and receive theirs into the halos (CommSend is
// non-blocking in both backends, so send-both-then-recv-both cannot deadlock). any(): OR each rank's
// "changed?" via a gather to rank 0 + broadcast (no collective primitive), which also barriers the round.
struct FlatComm {
  int r, ntiles; bool hasL, hasR; int tL, tR, tChg, tCont;
  template<class T> void exch(std::vector<T>& hL, std::vector<T>& hR,
                              const std::vector<T>& myL, const std::vector<T>& myR){
    if(hasR) c::CommSend(myR, r+1, tR);                 // my right edge -> right neighbour (its left halo)
    if(hasL) c::CommSend(myL, r-1, tL);                 // my left  edge -> left  neighbour (its right halo)
    if(hasL){ c::CommRecv(hL, r-1, tR); }               // left  neighbour's right edge (it sent with tR)
    if(hasR){ c::CommRecv(hR, r+1, tL); }               // right neighbour's left  edge (it sent with tL)
  }
  bool any(bool local){
    if(r==0){ bool a=local; uint8_t o;
      for(int t=1;t<ntiles;t++){ c::CommRecv(o,t,tChg); a = a||(o!=0); }
      uint8_t ao=a?1:0; for(int t=1;t<ntiles;t++) c::CommSend(ao,t,tCont); return a;
    } else { uint8_t lo=local?1:0; c::CommSend(lo,0,tChg); uint8_t ao; c::CommRecv(ao,0,tCont); return ao!=0; }
  }
};

int main(int argc, char **argv){
#ifdef DH_USE_MPI
  commt::CommInitMPI();                     // one rank per process; must precede any CommRank/Size
#endif
  if(argc!=4 && argc!=5){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level> <Split Cols (comma-sep)> [flat halo cap]\n";
    return -1;
  }
  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);
  const int         halo_cap    = (argc==5) ? std::stoi(argv[4]) : std::numeric_limits<int>::max();

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

#ifdef DH_USE_MPI
  if(ntiles != commt::CommSize()){
    if(commt::CommRank()==0)
      std::cerr<<"error: "<<ntiles<<" tiles from the split but "<<commt::CommSize()
               <<" MPI ranks; launch with `mpirun -np "<<ntiles<<"`\n";
    commt::CommFinalizeMPI();
    return 1;
  }
#endif

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
    rd::Array2D<dh_label_t> glab_pc;     // same, after conduit resolution (BOUNDARY cells resolved)
    rd::Array2D<int8_t>     gfix;        // tile flood flowdirs with the seam fix-up applied (pre-flats)
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
    dh::GetDepressionHierarchyPhaseAB<float,rd::Topology::D8>(dem, label, fd, deps, outlets, /*permit_without_baselevel_seed=*/true);
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

  // ORACLE conduit resolution (the stitch's boundary-graph pass, on the oracle tiles with a
  // full-grid drain): resolve every BOUNDARY cell to the real depression/ocean it drains into.
  // Phase 1 localwalk follows each tile's own fd to a local terminal or a seam exit; Phase 2
  // chase follows exits across tiles to a terminal. Produces the post-conduit global grid.
  rd::Array2D<dh_label_t> gLabel_oracle_pc = gLabel_oracle;
  {
    struct LW { bool terminal; dh_label_t label; int ex,ey; };
    const auto localwalk = [&](int t,int gx,int gy)->LW{
      const int x0=bounds[t];
      for(unsigned guard=0; guard<=oracle[t].label.size(); guard++){
        const dh_label_t ll = oracle[t].label(gx-x0,gy);
        if(ll!=BOUNDARY) return { true, gmap(ll,oracle[t].offset), 0,0 };
        const int8_t f = oracle[t].fd(gx-x0,gy);
        if(f!=rd::NO_FLOW){ gx+=rd::d8x[f]; gy+=rd::d8y[f]; }
        else { int nx,ny; bool toO;
          if(!drain_full(gx,gy,nx,ny,toO) || toO) return { true, OCEAN, 0,0 };
          return { false, OCEAN, nx, ny }; }
      }
      return { true, OCEAN, 0,0 };
    };
    std::map<std::pair<int,int>,LW> res;
    for(int t=0;t<ntiles;t++)
      for(int y=0;y<H;y++)
        for(int gx=bounds[t];gx<bounds[t+1];gx++)
          if(oracle[t].label(gx-bounds[t],y)==BOUNDARY) res[{gx,y}] = localwalk(t,gx,y);
    const auto chase = [&](int gx,int gy)->dh_label_t{
      for(unsigned guard=0; guard<=res.size()+1; guard++){
        const LW &r = res.at({gx,gy});
        if(r.terminal) return r.label;
        const int nt = tile_of(r.ex);
        const dh_label_t nl = oracle[nt].label(r.ex-bounds[nt], r.ey);
        if(nl!=BOUNDARY) return gmap(nl, oracle[nt].offset);
        gx=r.ex; gy=r.ey;
      }
      return OCEAN;
    };
    for(const auto &kv : res) gLabel_oracle_pc(kv.first.first, kv.first.second) = chase(kv.first.first, kv.first.second);
  }

  // ORACLE outlet database (the stitch's outlet step, on the oracle post-conduit grid): one
  // entry per depression pair, holding the outlet = the higher-elevation cell of the pair,
  // keeping the LOWEST such max per pair (first-seen wins ties). Intra-tile adjacencies, then
  // Barnes' HandleEdge across each seam. This is the object the distributed build must rebuild.
  using OutKey = std::pair<dh_label_t,dh_label_t>;
  using OutVal = std::pair<float,int64_t>;                 // (out_elev, out_cell)
  // Tie-break (WATCH note, resolved): on EQUAL outlet elevation keep the LOWER out_cell. This
  // reproduces the stitch's global-row-major "first-seen at lowest max" (first-seen in row-major
  // at a given elevation IS the lowest-index cell) as an order-INDEPENDENT rule, so the
  // distributed merge -- which visits pairs in rank order, not row-major -- picks the same cell.
  // It does not change serial/stitch output (row-major already yields the lowest index).
  const auto reduce = [](std::map<OutKey,OutVal> &db, dh_label_t la, dh_label_t lb,
                         float ea, int64_t ca, float eb, int64_t cb){
    if(la==lb) return;
    const float oe = (ea>=eb)? ea : eb; const int64_t oc = (ea>=eb)? ca : cb;
    const auto key = std::minmax(la,lb);
    const auto it = db.find({key.first,key.second});
    if(it==db.end() || oe < it->second.first || (oe==it->second.first && oc < it->second.second))
      db[{key.first,key.second}] = {oe,oc};
  };
  std::map<OutKey,OutVal> odb;
  // Skip a cell only if NoData AND not OCEAN -- a NoData-as-OCEAN cell is a real base level whose
  // adjacency to a land depression is a genuine basin->ocean outlet (mirrors the stitch fix 261bbcd;
  // the whole harness re-derives outlets from resolved labels, so the same NoData-skip lived here too).
  const auto skip_o = [&](int x,int y){ return full.isNoData(x,y) && gLabel_oracle_pc(x,y)!=OCEAN; };
  for(int y=0;y<H;y++)
    for(int x=0;x<W;x++){
      if(skip_o(x,y)) continue;
      for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
        if(!dx && !dy) continue;
        const int nx=x+dx, ny=y+dy;
        if(!full.inGrid(nx,ny) || skip_o(nx,ny)) continue;
        if(tile_of(nx)!=tile_of(x)) continue;             // cross-seam handled by HandleEdge
        reduce(odb, gLabel_oracle_pc(x,y), gLabel_oracle_pc(nx,ny),
               full(x,y), full.xyToI(x,y), full(nx,ny), full.xyToI(nx,ny));
      }
    }
  for(size_t b=1;b+1<bounds.size();b++){
    const int cA=bounds[b]-1, cB=bounds[b];
    for(int y=0;y<H;y++){
      if(skip_o(cA,y)) continue;
      for(int ny=y-1;ny<=y+1;ny++){
        if(ny<0 || ny>=H || skip_o(cB,ny)) continue;
        reduce(odb, gLabel_oracle_pc(cA,y), gLabel_oracle_pc(cB,ny),
               full(cA,y), full.xyToI(cA,y), full(cB,ny), full.xyToI(cB,ny));
      }
    }
  }

  // ---- DISTRIBUTED: one thread-rank per tile, BOUNDARY from a 1-column halo ----
  enum { TAG_L2R=1, TAG_R2L=2, TAG_COUNT=3, TAG_OFFSET=4,     // strip dirs; count->rank0; offset->rank r
         TAG_LW=5, TAG_EDGE=6, TAG_RESOLVED=7,                // conduit: Phase-1 recs / edge labels / results
         TAG_RES_LEFT=8, TAG_ODB_INTRA=9, TAG_ODB_EDGE=10,    // outlets: resolved edge strip / intra-db / edge-db
         TAG_DEPREC=11,                                       // depression records -> rank 0
         TAG_TREE_E=12, TAG_TREE_P=13,                        // Phase D: broadcast tree out_elev / parent
         TAG_MARG_C=14, TAG_MARG_E=15,                        // Phase D: reduce marginal cell_count / total_elev
         TAG_DISTSLICE=16,                                    // gather per-cell grid slices to rank 0 (MPI verify)
         TAG_FLAT_L=17, TAG_FLAT_R=18,                        // ENH-1 per-rank flats: seam field columns (L/R dir)
         TAG_FLAT_CHG=19, TAG_FLAT_CONT=20 };                 // ENH-1 per-rank flats: convergence gather / continue
  std::map<OutKey,OutVal>         outlet_db_dist;             // rank 0's merged outlet DB (read in verify)
  dh::DepressionHierarchy<float>  Gdist;                      // rank 0's assembled global hierarchy (leaves)
  dh_label_t                      n_global_r0 = 0;            // rank 0's global depression count (incl. ocean)
  auto rank_main = [&](){
    const int r  = c::CommRank();
    const int x0 = bounds[r], x1 = bounds[r+1], w = x1-x0;
    rd::Array2D<float>      dem   = extract_cols(full, x0, x1);   // this rank's own columns only
    rd::Array2D<dh_label_t> label = ocean_labels(dem, ocean_level);

    // Build this rank's edge strips and exchange with seam neighbours (non-blocking send
    // into the neighbour's inbox, then blocking recv -- no deadlock).
    const auto strip_of = [&](int lc){
      EdgeStrip s; s.elev.resize(H); s.ocean.resize(H); s.nodata.resize(H);
      for(int y=0;y<H;y++){ s.elev[y]=dem(lc,y); s.ocean[y]=(dem.isNoData(lc,y)||dem(lc,y)==ocean_level); s.nodata[y]=dem.isNoData(lc,y); }
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
    const auto drain_local = [&](int gx,int y,int &lx,int &ly,bool &to_ocean)->bool{
      const float focal = elev(gx,y);
      float best = std::numeric_limits<float>::infinity(); bool found=false; to_ocean=false;
      for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
        if(!dx && !dy) continue;
        const int nx=gx+dx, ny=y+dy;
        if(nx<0||nx>=W||ny<0||ny>=H) continue;
        if(is_ocean(nx,ny)){ to_ocean=true; continue; }
        const float e = elev(nx,ny);
        if(e >= focal) continue;
        if(e <= best){ best=e; lx=nx; ly=ny; found=true; }        // <= : higher-index tie wins
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
        int lx,ly; bool to_ocean;
        if(drain_local(gc,y,lx,ly,to_ocean) && !to_ocean && tile_of(lx)!=r){ label(lc,y)=BOUNDARY; nb++; }
      }

    rd::Array2D<int8_t> fd(dem.width(), dem.height(), rd::NO_FLOW);
    dh::DepressionHierarchy<float> deps; std::vector<dh::Outlet<float>> outlets;
    dh::GetDepressionHierarchyPhaseAB<float,rd::Topology::D8>(dem, label, fd, deps, outlets, /*permit_without_baselevel_seed=*/true);

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
      n_global_r0 = off;                                          // 1 + total depressions across tiles
    } else {
      c::CommSend(mycount, 0, TAG_COUNT);
      c::CommRecv(myoffset, 0, TAG_OFFSET);
    }

    // Map this rank's local labels into the global namespace (its own slice; no gather).
    rd::Array2D<dh_label_t> glab(w, H);
    for(int y=0;y<H;y++) for(int x=0;x<w;x++) glab(x,y) = gmap(label(x,y), myoffset);

    // ---- distributed conduit resolution (eng-doc section 2, v1 gather-and-resolve) ----
    // PHASE 1 (LOCAL): walk each BOUNDARY cell along this tile's own fd to a local terminal
    // (a global label) or, at a NO_FLOW seam seed, hand off across the seam (drain_local uses
    // the halo). Emits one O(boundary) LWRec per BOUNDARY cell. Also publish this tile's
    // edge-column global labels so rank 0 can read an EXIT target's label. Nothing but these
    // O(boundary) records leaves the rank.
    std::vector<LWRec>  mylw;
    std::vector<CellLab> myedge;
    const auto localwalk = [&](int gx,int gy)->LWRec{
      int cx=gx, cy=gy;
      for(unsigned guard=0; guard<=label.size(); guard++){
        const dh_label_t ll = label(cx-x0,cy);
        if(ll!=BOUNDARY) return { gx,gy, 1, gmap(ll,myoffset), 0,0 };   // local terminal
        const int8_t f = fd(cx-x0,cy);
        if(f!=rd::NO_FLOW){ cx+=rd::d8x[f]; cy+=rd::d8y[f]; }           // stays in this tile
        else { int nx,ny; bool toO;                                    // seam seed: cross the strip
          if(!drain_local(cx,cy,nx,ny,toO) || toO) return { gx,gy, 1, OCEAN, 0,0 };
          return { gx,gy, 0, 0, nx,ny }; }                             // EXIT to neighbour cell
      }
      return { gx,gy, 1, OCEAN, 0,0 };
    };
    for(int y=0;y<H;y++)                                               // ALL BOUNDARY cells: the label
      for(int gx=x0;gx<x1;gx++)                                        // spreads inward through PhaseAB,
        if(label(gx-x0,y)==BOUNDARY) mylw.push_back(localwalk(gx,y));  // not just the seeded seam cells
    for(int gc : seam_cols)                                            // edge-column labels (EXIT targets
      for(int y=0;y<H;y++) myedge.push_back({ gc, y, glab(gc-x0,y) }); // land on a neighbour's edge column)

    // PHASE 2 (rank 0): gather the O(boundary) records, chase each BOUNDARY cell to its
    // terminal, scatter the resolved labels back to their owning rank.
    std::vector<CellLab> myresolved;
    if(r==0){
      std::map<std::pair<int,int>,LWRec>     lwmap;
      std::map<std::pair<int,int>,dh_label_t> edgemap;
      const auto absorb = [&](const std::vector<LWRec> &lw, const std::vector<CellLab> &ed){
        for(const auto &q : lw) lwmap[{q.gx,q.gy}] = q;
        for(const auto &e : ed) edgemap[{e.gx,e.gy}] = e.label;
      };
      absorb(mylw, myedge);
      for(int t=1;t<ntiles;t++){
        std::vector<LWRec> lw; std::vector<CellLab> ed;
        c::CommRecv(lw, t, TAG_LW); c::CommRecv(ed, t, TAG_EDGE);
        absorb(lw, ed);
      }
      const auto chase = [&](int gx,int gy)->dh_label_t{
        int cx=gx, cy=gy;
        for(unsigned guard=0; guard<=lwmap.size()+1; guard++){
          const LWRec &q = lwmap.at({cx,cy});
          if(q.terminal) return q.label;
          const auto it = edgemap.find({q.ex,q.ey});                  // EXIT target's global label
          const dh_label_t nl = (it==edgemap.end()) ? OCEAN : it->second;
          if(nl!=BOUNDARY) return nl;                                 // terminal across the seam
          cx=q.ex; cy=q.ey;                                           // else chain on its LWRec
        }
        return OCEAN;
      };
      std::vector<std::vector<CellLab>> out(ntiles);                  // resolved labels per owning rank
      for(const auto &kv : lwmap){
        const int gx=kv.first.first, gy=kv.first.second;
        out[tile_of(gx)].push_back({ gx, gy, chase(gx,gy) });
      }
      for(int t=1;t<ntiles;t++) c::CommSend(out[t], t, TAG_RESOLVED);
      myresolved = std::move(out[0]);
    } else {
      c::CommSend(mylw,   0, TAG_LW);
      c::CommSend(myedge, 0, TAG_EDGE);
      c::CommRecv(myresolved, 0, TAG_RESOLVED);
    }

    // Overlay the resolved BOUNDARY labels onto this rank's global-label slice.
    rd::Array2D<dh_label_t> glab_pc = glab;
    for(const auto &c2 : myresolved) glab_pc(c2.gx-x0, c2.gy) = c2.label;

    // ---- distributed outlet set (eng-doc component 4) ----
    c::CommBarrier();                                            // separate the conduit gather from this one
    // Intra-tile outlet DB: adjacencies whose neighbour is in THIS tile (local: own columns only).
    std::map<OutKey,OutVal> myintra;
    // Skip a cell only if NoData AND not OCEAN (mirrors skip_o / the stitch fix 261bbcd): a
    // NoData-as-OCEAN cell's adjacency to a land depression is a real basin->ocean outlet.
    const auto skip_d = [&](int lx,int ly){ return dem.isNoData(lx,ly) && glab_pc(lx,ly)!=OCEAN; };
    for(int y=0;y<H;y++) for(int gx=x0;gx<x1;gx++){
      if(skip_d(gx-x0,y)) continue;
      for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
        if(!dx && !dy) continue;
        const int nx=gx+dx, ny=y+dy;
        if(nx<x0||nx>=x1||ny<0||ny>=H) continue;                // same-tile neighbours only
        if(skip_d(nx-x0,ny)) continue;
        reduce(myintra, glab_pc(gx-x0,y), glab_pc(nx-x0,ny),
               dem(gx-x0,y), (int64_t)y*W+gx, dem(nx-x0,ny), (int64_t)ny*W+nx);
      }
    }
    // HandleEdge across seams: send my LEFT edge column (resolved) to r-1; each rank runs the
    // match for its RIGHT seam using the neighbour's LEFT-edge strip. O(H) per seam.
    if(r>0){
      ResStrip s; s.label.resize(H); s.elev.resize(H); s.nodata.resize(H);
      for(int y=0;y<H;y++){ s.label[y]=glab_pc(0,y); s.elev[y]=dem(0,y); s.nodata[y]=dem.isNoData(0,y); }
      c::CommSend(s, r-1, TAG_RES_LEFT);
    }
    std::map<OutKey,OutVal> myedgedb;
    if(r<ntiles-1){
      ResStrip nbr; c::CommRecv(nbr, r+1, TAG_RES_LEFT);
      const int cA=x1-1, cB=bounds[r+1];                        // own right edge / neighbour left edge (x1)
      for(int y=0;y<H;y++){
        if(skip_d(cA-x0,y)) continue;
        for(int ny=y-1;ny<=y+1;ny++){
          if(ny<0||ny>=H) continue;
          if(nbr.nodata[ny] && nbr.label[ny]!=OCEAN) continue;  // NoData neighbour kept only if OCEAN
          reduce(myedgedb, glab_pc(cA-x0,y), nbr.label[ny],
                 dem(cA-x0,y), (int64_t)y*W+cA, nbr.elev[ny], (int64_t)ny*W+cB);
        }
      }
    }
    // Gather per-rank DBs to rank 0; merge all intra (rank order), then all edge (rank order),
    // matching the stitch's "all intra-tile, then all seams" outlet order.
    const auto to_vec=[&](const std::map<OutKey,OutVal> &db){
      std::vector<ORec> v; v.reserve(db.size());
      for(const auto &kv : db) v.push_back({ kv.first.first, kv.first.second, kv.second.first, kv.second.second });
      return v;
    };
    std::vector<ORec> intra_vec=to_vec(myintra), edge_vec=to_vec(myedgedb);
    if(r==0){
      std::vector<std::vector<ORec>> intra_all(ntiles), edge_all(ntiles);
      intra_all[0]=intra_vec; edge_all[0]=edge_vec;
      for(int t=1;t<ntiles;t++){ c::CommRecv(intra_all[t],t,TAG_ODB_INTRA); c::CommRecv(edge_all[t],t,TAG_ODB_EDGE); }
      std::map<OutKey,OutVal> gdb;
      const auto mrg=[&](const ORec &rc){                       // same tie-break as reduce(): lower cell wins
        const auto it=gdb.find({rc.depa,rc.depb});
        if(it==gdb.end() || rc.oelev < it->second.first
                         || (rc.oelev==it->second.first && rc.ocell < it->second.second))
          gdb[{rc.depa,rc.depb}]={rc.oelev,rc.ocell};
      };
      for(int t=0;t<ntiles;t++) for(const auto &rc:intra_all[t]) mrg(rc);
      for(int t=0;t<ntiles;t++) for(const auto &rc:edge_all[t]) mrg(rc);
      outlet_db_dist = std::move(gdb);
    } else {
      c::CommSend(intra_vec, 0, TAG_ODB_INTRA);
      c::CommSend(edge_vec,  0, TAG_ODB_EDGE);
    }

    // ---- gather depression records to rank 0 (eng-doc component 7) ----
    // Each rank ships its leaf depressions with dep_label and pit_cell remapped to global.
    // These O(#local depressions) records are the tree's raw material -- the one inherently
    // global object (the meta-tree has no per-tile home), centralized on rank 0 for Phase C.
    std::vector<dh::Depression<float>> myrecs;
    for(dh_label_t k=1;k<deps.size();k++){
      auto d = deps[k];
      d.dep_label = gmap(k, myoffset);
      if(d.pit_cell != dh::NO_VALUE){
        int lx,ly; dem.iToxy(d.pit_cell, lx, ly);
        d.pit_cell = (dh::flat_c_idx)ly*W + (x0+lx);              // tile-local -> global row-major index
      }
      myrecs.push_back(d);
    }
    if(r==0){
      Gdist = dh::DepressionHierarchy<float>(n_global_r0);
      Gdist[0] = deps[0]; Gdist[0].dep_label = 0;                 // the shared ocean node
      const auto place=[&](const std::vector<dh::Depression<float>> &v){ for(const auto &d:v) Gdist[d.dep_label]=d; };
      place(myrecs);
      for(int t=1;t<ntiles;t++){ std::vector<dh::Depression<float>> v; c::CommRecv(v,t,TAG_DEPREC); place(v); }
    } else {
      c::CommSend(myrecs, 0, TAG_DEPREC);
    }

    // ---- Phase C (rank 0, grid-free assembly) + distributed Phase D (marginal volumes) ----
    // Rank 0 assembles the hierarchy from the gathered records + distributed outlet set, then
    // broadcasts the light tree (out_elev + parent per node). Each rank computes the marginal
    // (cell_count, total_elevation) contribution of ITS OWN cells -- the exact per-cell walk-up
    // of CalculateMarginalVolumes, but tile-local -- and the O(#deps) partials reduce to rank 0,
    // which applies them and runs the grid-free CalculateTotalVolumes. No rank but the owner
    // reads any cell, so the last full-grid dependency (Phase D) is now distributed.
    c::CommBarrier();
    std::vector<float>      tree_out_elev;
    std::vector<dh_label_t> tree_parent;
    if(r==0){
      std::vector<dh::Outlet<float>> outlets;
      for(const auto &kv : outlet_db_dist){
        dh::Outlet<float> o; o.depa=kv.first.first; o.depb=kv.first.second;
        o.out_elev=kv.second.first; o.out_cell=kv.second.second; outlets.push_back(o);
      }
      dh::GetDepressionHierarchyPhaseC<float>(Gdist, outlets);   // grid-free; grows Gdist with meta nodes
      const dh_label_t T = Gdist.size();
      tree_out_elev.resize(T); tree_parent.resize(T);
      for(dh_label_t i=0;i<T;i++){ tree_out_elev[i]=Gdist[i].out_elev; tree_parent[i]=Gdist[i].parent; }
      for(int t=1;t<ntiles;t++){ c::CommSend(tree_out_elev,t,TAG_TREE_E); c::CommSend(tree_parent,t,TAG_TREE_P); }
    } else {
      c::CommRecv(tree_out_elev, 0, TAG_TREE_E);
      c::CommRecv(tree_parent,   0, TAG_TREE_P);
    }
    const dh_label_t T = tree_out_elev.size();
    std::vector<uint32_t> dcount(T, 0);
    std::vector<double>   delev (T, 0.0);
    for(int y=0;y<H;y++) for(int gx=x0;gx<x1;gx++){             // this rank's own cells only
      if(dem.isNoData(gx-x0,y)) continue;
      const float me = dem(gx-x0,y);
      dh_label_t cl = glab_pc(gx-x0,y);
      while(cl!=OCEAN && me > tree_out_elev[cl]) cl = tree_parent[cl];   // walk up to the depression it belongs to
      if(cl==OCEAN) continue;
      dcount[cl]++; delev[cl]+=me;
    }
    if(r==0){
      for(int t=1;t<ntiles;t++){
        std::vector<uint32_t> dc; std::vector<double> de;
        c::CommRecv(dc,t,TAG_MARG_C); c::CommRecv(de,t,TAG_MARG_E);
        for(dh_label_t i=0;i<T;i++){ dcount[i]+=dc[i]; delev[i]+=de[i]; }
      }
      for(dh_label_t i=0;i<T;i++){ Gdist[i].cell_count += dcount[i]; Gdist[i].total_elevation += delev[i]; }
      dh::CalculateTotalVolumes<float>(Gdist);                 // grid-free rollup
    } else {
      c::CommSend(dcount,0,TAG_MARG_C); c::CommSend(delev,0,TAG_MARG_E);
    }

    // ---- per-cell flowdir seam fix-up (eng-doc section 6, non-flat crossings) ----
    // A tile flood cannot point a cell across its own boundary, so a cell whose true drainage
    // exits the tile is left NO_FLOW (a land seam seed) or points at the wrong ocean cell. Every
    // other cell already matches serial. Restore serial's choice using ONLY the 1-column halo:
    // a land seam seed points to its lowest cross-tile neighbour at/below its pit (highest-index
    // tie); a sea-draining cell points to its highest-index adjacent ocean cell. O(boundary).
    const auto isnodata=[&](int gx,int y)->bool{
      if(gx>=x0 && gx<x1) return dem.isNoData(gx-x0,y);
      if(gx==x0-1)        return haloL.nodata[y];
      return haloR.nodata[y];
    };
    const auto gidx=[&](int gx,int y)->int64_t{ return (int64_t)y*W+gx; };
    const auto dir_to=[&](int dx,int dy)->int8_t{ for(int n=1;n<=8;n++) if(rd::d8x[n]==dx && rd::d8y[n]==dy) return (int8_t)n; return rd::NO_FLOW; };
    rd::Array2D<int8_t> gfix = fd;                               // start from the tile's flood flowdirs
    for(int y=0;y<H;y++) for(int lx=0;lx<w;lx++){
      const int gx=x0+lx;
      if(fd(lx,y)==rd::NO_FLOW && !dem.isNoData(lx,y) && label(lx,y)!=OCEAN){
        const float fe=dem(lx,y); const int64_t fi=gidx(gx,y);
        int bnx=0,bny=0; int64_t bi=-1; float be=std::numeric_limits<float>::infinity();
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
          if(!dx && !dy) continue;
          const int nx=gx+dx, ny=y+dy;
          if(nx<0||nx>=W||ny<0||ny>=H || isnodata(nx,ny)) continue;
          if(!(nx<x0 || nx>=x1)) continue;                       // cross-seam neighbours only
          const float e=elev(nx,ny); if(e>fe) continue;          // at or below the pit
          const int64_t idx=gidx(nx,ny);
          if(e<be || (e==be && idx>bi)){ be=e; bi=idx; bnx=nx; bny=ny; }
        }
        if(bi>=0 && (be<fe || bi>fi)) gfix(lx,y)=dir_to(bnx-gx, bny-y);
      } else if(label(lx,y)==OCEAN && !dem.isNoData(lx,y)){       // sea-draining cell
        int64_t bi=-1; int bnx=0,bny=0;
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
          if(!dx && !dy) continue;
          const int nx=gx+dx, ny=y+dy;
          if(nx<0||nx>=W||ny<0||ny>=H || !is_ocean(nx,ny)) continue;
          const int64_t idx=gidx(nx,ny);
          if(idx>bi){ bi=idx; bnx=nx; bny=ny; }
        }
        if(bi>=0 && (bnx<x0 || bnx>=x1)) gfix(lx,y)=dir_to(bnx-gx, bny-y);
      }
    }

    // ENH-1: resolve THIS tile's flats per-rank (bit-identical to serial), holding only the tile + a
    // 1-column halo -- three monotone relaxations with a per-round seam exchange + convergence all-reduce.
    { FlatComm fc{r, ntiles, r>0, r<ntiles-1, TAG_FLAT_L, TAG_FLAT_R, TAG_FLAT_CHG, TAG_FLAT_CONT};
      resolve_flat_flowdirs_rank(dem, gfix, fc); }

    dist[r].label = std::move(label);                            // ranks write disjoint slots
    dist[r].gfix  = std::move(gfix);
    dist[r].fd    = std::move(fd);
    dist[r].glab  = std::move(glab);
    dist[r].glab_pc = std::move(glab_pc);
    dist[r].nboundary = nb;
    dist[r].ndep  = mycount;
    dist[r].offset= myoffset;
    c::CommBarrier();
  };
  c::CommInit(ntiles, rank_main);

#ifdef DH_USE_MPI
  // Real MPI: each rank has only its own dist[] slot (separate address spaces). Gather every rank's
  // per-cell grid slices to rank 0 so the FULL per-cell verification below runs under real MPI, not
  // just the tree. Non-zero ranks ship their slice and exit; rank 0 fills dist[] and verifies.
  if(commt::CommRank()!=0){
    const auto &R = dist[commt::CommRank()];
    DistSlice s; s.w=R.label.width(); s.h=R.label.height(); const int nn=s.w*s.h;
    s.nboundary=R.nboundary; s.ndep=R.ndep; s.offset=R.offset;
    s.label.resize(nn); s.glab.resize(nn); s.glab_pc.resize(nn); s.fd.resize(nn); s.gfix.resize(nn);
    for(int i=0;i<nn;i++){ s.label[i]=R.label(i); s.glab[i]=R.glab(i); s.glab_pc[i]=R.glab_pc(i); s.fd[i]=R.fd(i); s.gfix[i]=R.gfix(i); }
    c::CommSend(s, 0, TAG_DISTSLICE);
    commt::CommFinalizeMPI();
    return 0;
  }
  for(int t=1;t<ntiles;t++){
    DistSlice s; c::CommRecv(s, t, TAG_DISTSLICE);
    auto &R = dist[t]; const int nn=s.w*s.h;
    R.nboundary=s.nboundary; R.ndep=s.ndep; R.offset=s.offset;
    R.label  =rd::Array2D<dh_label_t>(s.w,s.h); R.glab=rd::Array2D<dh_label_t>(s.w,s.h);
    R.glab_pc=rd::Array2D<dh_label_t>(s.w,s.h); R.fd  =rd::Array2D<int8_t>(s.w,s.h);
    R.gfix   =rd::Array2D<int8_t>(s.w,s.h);
    for(int i=0;i<nn;i++){ R.label(i)=s.label[i]; R.glab(i)=s.glab[i]; R.glab_pc(i)=s.glab_pc[i]; R.fd(i)=s.fd[i]; R.gfix(i)=s.gfix[i]; }
  }
#endif

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

  // ---- conduit check: the resolved global label grid (BOUNDARY cells now resolved) must
  // equal the stitch's post-conduit grid, cell for cell. ----
  long cdiff=0, nboundary=0;
  rd::Array2D<dh_label_t> gLabel_dist_pc(W, H, OCEAN);
  for(int t=0;t<ntiles;t++){
    const int x0=bounds[t], x1=bounds[t+1];
    for(int y=0;y<H;y++) for(int x=x0;x<x1;x++)
      gLabel_dist_pc(x,y) = dist[t].glab_pc(x-x0,y);
  }
  for(unsigned i=0;i<gLabel_oracle_pc.size();i++){
    if(gLabel_oracle(i)==BOUNDARY) nboundary++;          // cells the conduit had to resolve
    if(gLabel_oracle_pc(i)!=gLabel_dist_pc(i)) cdiff++;
  }
  const bool conduit_ok = (cdiff==0);
  std::cout<<(conduit_ok ? "MPI-CONDUIT-MATCH " : "MPI-CONDUIT-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" conduit_diff="<<cdiff<<"/"<<cells
           <<" boundary_resolved="<<nboundary<<"\n";

  // ---- outlet check: the merged distributed outlet DB must equal the stitch's, pair for pair.
  // Break out cell-only diffs (same pair, same outlet elevation, different equal-elev cell) --
  // those are order-dependent tie choices that may or may not perturb the final tree. ----
  long pair_missing=0, elev_diff=0, cell_only=0, extra=0;
  for(const auto &kv : odb){
    const auto it = outlet_db_dist.find(kv.first);
    if(it==outlet_db_dist.end()){ pair_missing++; continue; }
    if(it->second.first != kv.second.first)       elev_diff++;
    else if(it->second.second != kv.second.second) cell_only++;
  }
  for(const auto &kv : outlet_db_dist) if(!odb.count(kv.first)) extra++;
  const bool outlet_ok = (pair_missing==0 && elev_diff==0 && cell_only==0 && extra==0);
  std::cout<<(outlet_ok ? "MPI-OUTLET-MATCH " : "MPI-OUTLET-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" pairs="<<odb.size()
           <<" missing="<<pair_missing<<" extra="<<extra
           <<" elev_diff="<<elev_diff<<" cell_only_diff="<<cell_only<<"\n";

  // ---- tree check: rank_main already assembled the hierarchy (grid-free Phase C) and computed
  // volumes with the DISTRIBUTED Phase D (each rank's own cells, reduced to rank 0) -- no rank
  // read a foreign tile's interior. Here we only contract seam artifacts and require the
  // canonical signature (which includes per-node cell_count and dep_vol) to equal serial's. ----
  const int n_collapsed = CollapseSeamArtifacts(Gdist, full, bounds);

  // serial ground truth (tree + per-cell flowdirs)
  auto s_label = ocean_labels(full, ocean_level);
  rd::Array2D<int8_t> s_fd(W, H, rd::NO_FLOW);
  auto S = dh::GetDepressionHierarchy<float,rd::Topology::D8>(full, s_label, s_fd);
  resolve_flat_flowdirs(full, s_fd);                     // serial uses the same deterministic flat routing

  const std::string sig_dist = dhtest::canonicalize(Gdist);
  const std::string sig_serial = dhtest::canonicalize(S);
  const auto iv_d = dhtest::invariants(Gdist);
  const auto iv_s = dhtest::invariants(S);
  const bool tree_ok = (sig_dist==sig_serial);
  std::cout<<(tree_ok ? "MPI-TREE-MATCH " : "MPI-TREE-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" collapsed="<<n_collapsed
           <<" nodes(serial="<<iv_s.n_nodes<<" dist="<<iv_d.n_nodes<<")"
           <<" total_dep_vol(serial="<<iv_s.total_dep_vol<<" dist="<<iv_d.total_dep_vol<<")\n";

  // ---- decomposition-correctness diagnostic: the NODE COUNT ----
  // The depression tree is split-invariant: a correct build has the SAME set of depressions no matter
  // where the tile seams fall, so its node count is independent of the decomposition. An unequal node
  // count is therefore a NECESSARY-condition symptom that the split changed the tree -- weaker than the
  // full canonical signature but a real structural signal. It flags TWO things, which the diff must
  // separate by hand:
  //   (1) genuine seam artifacts the collapse pass failed to contract -- a spurious extra depression
  //       (Corsica coastal leaf pre-NoData-fix; testdem8 rim fragment pre-Pass-B2; the residual
  //       split-10 kerry cases). These are real bugs to chase; the collapse should make node count agree.
  //   (2) the STRUCTURAL sub-class of the PhaseCD tie-break: at a TIED outlet elevation two basins may be
  //       rebuilt either as co-equal children of a META or with one ocean_linked into the other -- SAME
  //       depressions and volume, but meta-vs-ocean_linked changes the node count (kerry_test2: serial
  //       meta+2 leaves vs stitch 1 leaf + 1 ocean_linked, 4 vs 3 nodes). This CHANGES serial output ->
  //       Richard-coordinated, acceptable under the volume-correct+valid-tree bar.
  // (So node count does NOT cleanly isolate "decomposition bug" from all tie-break noise: only the
  // pure-ORDERING tie-break sub-class preserves node count -- those show DIFFER but DECOMP-CORRECT, e.g.
  // kerry_test4 splits 3/8. The meta-vs-ocean_linked sub-class does not.) At GEBCO scale there is no
  // serial to diff; the same invariant then powers a serial-free self-check -- build at two splits, node
  // counts must agree (a mismatch there is class (1) or (2), never a decomposition that is actually right).
  const bool decomp_ok = (iv_d.n_nodes == iv_s.n_nodes);
  std::cout<<(decomp_ok ? "MPI-DECOMP-CORRECT " : "MPI-DECOMP-INCORRECT ")<<in_name
           <<" ranks="<<ntiles<<" nodes(serial="<<iv_s.n_nodes<<" dist="<<iv_d.n_nodes<<")\n";

  // Volume-correctness verdict (VOL-MATCH/VOL-DIFFER): total depression volume vs serial -- the acceptance
  // bar for cases NOT bit-identical to serial (bowl-interior ENH-2 splits, the tie-break class): they are
  // volume-correct + valid-tree by design even when the tree DIFFERs. Small relative tolerance absorbs the
  // distributed PhaseD's per-rank summation order; an open depression (inf) or a lost basin -> DIFFER.
  const double v_s = iv_s.total_dep_vol, v_d = iv_d.total_dep_vol;
  const bool vol_ok = std::isfinite(v_s) && std::isfinite(v_d) &&
                      std::fabs(v_d - v_s) <= 1e-6 * std::max(1.0, std::fabs(v_s));
  std::cout<<(vol_ok ? "MPI-VOL-MATCH " : "MPI-VOL-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" total_dep_vol(serial="<<v_s<<" dist="<<v_d<<")\n";

  // ---- flowdir check: each rank now resolves ITS flats per-rank (ENH-1, resolve_flat_flowdirs_rank in
  // rank_main above) using a 1-column seam exchange + convergence all-reduce -- so no rank ever holds a
  // full-grid field. Here we only ASSEMBLE the per-rank results for the differential check; the flat pass
  // no longer runs centrally. Must equal serial (resolve_flat_flowdirs on the full grid), cell for cell. ----
  rd::Array2D<int8_t> gFix(W, H, rd::NO_FLOW);
  for(int t=0;t<ntiles;t++){
    const int x0=bounds[t], x1=bounds[t+1];
    for(int y=0;y<H;y++) for(int x=x0;x<x1;x++)
      gFix(x,y) = dist[t].gfix(x-x0,y);
  }
  (void)halo_cap;                              // option 2 (now per-rank) supersedes option 3's cap; arg kept for compat
  long fd_land=0, fd_diff=0;
  rd::Array2D<int8_t> sdf(W,H,rd::NO_FLOW); rd::d8_flow_directions(full, sdf);   // to split diffs flat vs non-flat
  long fd_flat=0;                             // diffs at FLAT cells = the per-rank flat resolution's own error
  for(unsigned i=0;i<full.size();i++){
    if(full.isNoData(i)) continue;
    fd_land++;
    if(gFix(i)!=s_fd(i)){ fd_diff++; if(sdf(i)==rd::NO_FLOW) fd_flat++; }
  }
  const bool flowdir_ok = (fd_diff==0);
  std::cout<<(flowdir_ok ? "MPI-FLOWDIR-MATCH " : "MPI-FLOWDIR-DIFFER ")<<in_name
           <<" ranks="<<ntiles<<" fd_diff="<<fd_diff<<"/"<<fd_land<<" flat_diff="<<fd_flat<<"\n";

  const int rc = (phaseab_ok && remap_ok && conduit_ok && outlet_ok && tree_ok && flowdir_ok) ? 0 : 1;
#ifdef DH_USE_MPI
  commt::CommFinalizeMPI();
#endif
  return rc;
}
