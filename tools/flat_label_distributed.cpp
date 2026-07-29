// Lever C core (in-process de-risk): reproduce serial's flat-label PARTITION in the DISTRIBUTED shape and
// validate the ASSEMBLED global partition equals serial across tilings. No MPI yet; proves the algorithm
// before the plumbing (dephier_mpi). Exits 0 iff DIST-PARTITION-MATCH.
//
//   per tile t owning columns [b_t, b_{t+1}): replay the flood label partition over an adaptive window
//     [b_t-halo, b_{t+1}+halo), doubling halo until the OWNED columns' pit-labels are stable for TWO
//     consecutive doublings (one stable step can be a coincidence on a wide flat -- observed). Each tile
//     reads only its columns + halo -> O(N/P)+O(halo.boundary). A giant flat spanning many tiles drives the
//     halo up (toward the full grid for an EXACT result); a real build caps it (option 3) and accepts a
//     valid-but-maybe-not-identical partition beyond the cap, exactly as the flowdir halo does.
//   label = the basin's PIT GLOBAL CELL INDEX. This is namespace-free and consistent across tiles, so a
//     straddling basin carries the SAME label from every tile whose halo reaches its pit -- NO seam stitch /
//     union-find is needed; the assembled partition is just each cell's owning-tile pit-label.
//   compare: assembled vs serial's partition (namespace-free bijection).
#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;
using dh::dh_label_t; using dh::OCEAN; using dh::NO_DEP;

int main(int argc, char**argv){
  if(argc<4){ std::cerr<<"usage: "<<argv[0]<<" <dem> <ocean_level> <split_cols csv>\n"; return 2; }
  rd::Array2D<float> dem(argv[1]);
  const float ocean_level = std::stod(argv[2]);
  const int W=dem.width(), H=dem.height();
  const auto is_ocean=[&](int x,int y){ return dem.isNoData(x,y) || dem(x,y)==ocean_level; };

  std::vector<int> bounds={0};
  { std::string s=argv[3]; size_t p=0; while(true){ size_t q=s.find(',',p); bounds.push_back(std::stoi(s.substr(p,q-p))); if(q==std::string::npos)break; p=q+1; } }
  bounds.push_back(W);
  const int ntiles=bounds.size()-1;
  const auto tile_of=[&](int x){ int t=0; while(t+1<(int)bounds.size()&&x>=bounds[t+1])t++; return t; };

  // Windowed pit-up flood-label replay over columns [xlo,xhi), GLOBAL-index pop order. Fresh namespace.
  const auto replay_window=[&](int xlo,int xhi)->rd::Array2D<dh_label_t>{
    xlo=std::max(0,xlo); xhi=std::min(W,xhi);
    rd::Array2D<dh_label_t> r(W,H,NO_DEP);
    std::map<float,std::vector<long>> buckets;
    const auto inwin=[&](int x){ return x>=xlo&&x<xhi; };
    const auto push=[&](int x,int y){ buckets[dem(x,y)].push_back((long)y*W+x); };
    for(int y=0;y<H;y++) for(int x=xlo;x<xhi;x++){
      if(is_ocean(x,y)){ r(x,y)=OCEAN; push(x,y); continue; }
      const float e=dem(x,y); bool lower=false;
      for(int dy=-1;dy<=1&&!lower;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=x+dx,ny=y+dy;
        if(inwin(nx)&&dem.inGrid(nx,ny)&&dem(nx,ny)<e){ lower=true; break; } }
      if(!lower) push(x,y);
    }
    while(!buckets.empty()){
      auto it=buckets.begin(); const float e=it->first;
      std::vector<long> cur=std::move(it->second); buckets.erase(it);
      std::sort(cur.begin(),cur.end());
      while(!cur.empty()){ const long ci=cur.back(); cur.pop_back(); const int cx=ci%W, cy=ci/W;
        // Label a basin by its PIT's GLOBAL cell index -- namespace-free and consistent across tiles, so a
        // straddling basin gets the SAME label from every tile that reaches its pit (no union-find needed).
        dh_label_t cl=r(cx,cy); if(cl==NO_DEP){ cl=(dh_label_t)ci; r(cx,cy)=cl; }
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=cx+dx,ny=cy+dy;
          if(!inwin(nx)||!dem.inGrid(nx,ny)||r(nx,ny)!=NO_DEP) continue;
          r(nx,ny)=cl; if(dem(nx,ny)==e) cur.push_back((long)ny*W+nx); else buckets[dem(nx,ny)].push_back((long)ny*W+nx); }
      }
    }
    return r;
  };

  // Per-tile adaptive-halo converged window labelling.
  std::vector<rd::Array2D<dh_label_t>> Wt(ntiles);
  std::vector<int> used_halo(ntiles,0);
  // owned pit-labels unchanged between two window labellings, over columns [o0,o1)? (pit index is GLOBAL,
  // so a raw equality test is namespace-free.)
  const auto owned_stable=[&](const rd::Array2D<dh_label_t>&P,const rd::Array2D<dh_label_t>&C,int o0,int o1){
    for(int y=0;y<H;y++) for(int x=o0;x<o1;x++){ if(is_ocean(x,y))continue; if(P(x,y)!=C(x,y)) return false; }
    return true; };
  const int hmax=W;
  const bool full   = std::getenv("DH_FULL")!=nullptr;         // diag: window = full grid per tile (exact)
  const int  need   = std::getenv("DH_STABLE1")?1:2;           // consecutive stable doublings required (2 default)
  for(int t=0;t<ntiles;t++){
    const int o0=bounds[t], o1=bounds[t+1];
    if(full){ Wt[t]=replay_window(0,W); used_halo[t]=W; continue; }
    int h=std::min(2,hmax);
    rd::Array2D<dh_label_t> prev = replay_window(o0-h,o1+h);   // start at halo 2 (flowdir pattern)
    int stable_run=0;
    while(h<hmax){
      const int h2=std::min(h*2,hmax);
      rd::Array2D<dh_label_t> cur = replay_window(o0-h2,o1+h2);
      const bool same = owned_stable(prev,cur,o0,o1);
      prev=std::move(cur); h=h2;                               // ALWAYS advance to the larger window
      stable_run = same ? stable_run+1 : 0;
      if(stable_run>=need) break;
    }
    used_halo[t]=h;
    Wt[t]=std::move(prev);
  }

  // Assembled global partition = each cell's OWNING tile's pit-index (no stitch needed: the pit index is a
  // GLOBAL cell, so straddling basins carry the same label from every tile). Compare to serial (bijection).
  const auto plab=[&](int x,int y){ return Wt[tile_of(x)](x,y); };

  rd::Array2D<dh_label_t> s_label(W,H,NO_DEP);
  for(int y=0;y<H;y++) for(int x=0;x<W;x++) if(is_ocean(x,y)) s_label(x,y)=OCEAN;
  rd::Array2D<int8_t> s_fd(W,H,rd::NO_FLOW);
  dh::GetDepressionHierarchy<float,rd::Topology::D8>(dem, s_label, s_fd);

  std::map<dh_label_t,dh_label_t> s2r, r2s; long land=0, conf=0;
  for(int y=0;y<H;y++) for(int x=0;x<W;x++){ if(is_ocean(x,y))continue; land++;
    const dh_label_t s=s_label(x,y), r=plab(x,y);
    auto ia=s2r.find(s); if(ia==s2r.end())s2r[s]=r; else if(ia->second!=r)conf++;
    auto ib=r2s.find(r); if(ib==r2s.end())r2s[r]=s; else if(ib->second!=s)conf++; }
  int maxhalo=0; for(int h:used_halo) maxhalo=std::max(maxhalo,h);
  std::cout<<argv[1]<<" splits="<<argv[3]<<": "<<(conf==0?"DIST-PARTITION-MATCH":"DIST-PARTITION-DIFFER")
           <<" (land="<<land<<" conflicts="<<conf<<" max_halo="<<maxhalo<<" ntiles="<<ntiles<<")\n";
  return conf==0?0:1;
}
