// FEASIBILITY / regression guard (lever C, SPLIT_INVARIANT_FLATS_PLAN.md): can serial's flat-label PARTITION
// be reproduced across a seam by a BOUNDED local replay, and what halo does it need? Mirrors the proven
// flowdir adaptive-halo approach, for LABELS. Exits 0 iff every seam reproduces serial's owned partition with
// a halo strictly smaller than the full grid (bounded locality) -- the property lever C's dephier_mpi
// implementation relies on. The conclusion (form 2 = adaptive-halo ordered replay, O(flat-extent) halo, NOT
// clean O(boundary) relaxation because the partition is order-dependent) is recorded in the plan.
//
// Method: `replay_window(xlo,xhi)` runs the full pit-up flood-label replay RESTRICTED to columns [xlo,xhi),
// using GLOBAL row-major indices for the (elevation, index) pop order (so order matches the whole-grid flood
// within a column window) and treating the window edges as grid edges. The reference is the whole-grid
// window [0,W) (== serial's partition, proven). For a seam at column s, the LEFT tile owns [0,s): grow a halo
// H so the tile replays [0, s+H); we find the minimal H for which the OWNED cells' partition matches the
// reference. Symmetric for the RIGHT tile [s,W) replaying [s-H, W). Reports the halo per seam and correlates
// it with the widest seam-straddling connected flat. If a bounded H suffices, the halo-cap bounded-replay
// form is viable (like flowdir option-3); a required H that scales with flat width is that same known
// tradeoff, not a blocker.
#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;
using dh::dh_label_t; using dh::OCEAN; using dh::NO_DEP;

int main(int argc, char**argv){
  if(argc<4){ std::cerr<<"usage: "<<argv[0]<<" <dem> <ocean_level> <seam> [seam2 ...]\n"; return 2; }
  rd::Array2D<float> dem(argv[1]);
  const float ocean_level = std::stod(argv[2]);
  const int W=dem.width(), H=dem.height();
  const auto is_ocean=[&](int x,int y){ return dem.isNoData(x,y) || dem(x,y)==ocean_level; };

  // Windowed pit-up flood-label replay over columns [xlo,xhi). GLOBAL indices order the buckets.
  const auto replay_window=[&](int xlo,int xhi)->rd::Array2D<dh_label_t>{
    rd::Array2D<dh_label_t> r(W,H,NO_DEP);
    std::map<float,std::vector<long>> buckets;
    const auto inwin=[&](int x){ return x>=xlo && x<xhi; };
    const auto push=[&](int x,int y){ buckets[dem(x,y)].push_back((long)y*W+x); };
    for(int y=0;y<H;y++) for(int x=xlo;x<xhi;x++){
      if(is_ocean(x,y)){ r(x,y)=OCEAN; push(x,y); continue; }
      const float e=dem(x,y); bool lower=false;
      for(int dy=-1;dy<=1&&!lower;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=x+dx,ny=y+dy;
        if(inwin(nx)&&dem.inGrid(nx,ny)&&dem(nx,ny)<e){ lower=true; break; } }
      if(!lower) push(x,y);
    }
    dh_label_t next=1;
    while(!buckets.empty()){
      auto it=buckets.begin(); const float e=it->first;
      std::vector<long> cur=std::move(it->second); buckets.erase(it);
      std::sort(cur.begin(),cur.end());                      // global index asc; pop_back = highest first
      while(!cur.empty()){
        const long ci=cur.back(); cur.pop_back(); const int cx=ci%W, cy=ci/W;
        dh_label_t cl=r(cx,cy);
        if(cl==NO_DEP){ cl=next++; r(cx,cy)=cl; }
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=cx+dx,ny=cy+dy;
          if(!inwin(nx)||!dem.inGrid(nx,ny)||r(nx,ny)!=NO_DEP) continue;
          r(nx,ny)=cl; if(dem(nx,ny)==e) cur.push_back((long)ny*W+nx); else buckets[dem(nx,ny)].push_back((long)ny*W+nx); }
      }
    }
    return r;
  };

  // Partition agreement (namespace-free bijection) over a column range [ox0,ox1) of owned land cells.
  const auto owned_matches=[&](const rd::Array2D<dh_label_t>&ref,const rd::Array2D<dh_label_t>&L,int ox0,int ox1)->long{
    std::map<dh_label_t,dh_label_t> a,b; long conf=0;
    for(int y=0;y<H;y++) for(int x=ox0;x<ox1;x++){
      if(is_ocean(x,y)) continue;
      const dh_label_t s=ref(x,y), r=L(x,y);
      auto ia=a.find(s); if(ia==a.end()) a[s]=r; else if(ia->second!=r) conf++;
      auto ib=b.find(r); if(ib==b.end()) b[r]=s; else if(ib->second!=s) conf++;
    }
    return conf;
  };

  const auto ref = replay_window(0,W);   // == serial's partition (proven)

  int seams=0, bounded=0;
  for(int ai=3; ai<argc; ai++){
    const int s=std::stoi(argv[ai]);
    if(s<=0||s>=W){ std::cout<<"seam "<<s<<" out of range\n"; continue; }
    // widest seam-straddling connected flat (same-elevation component touching both sides of column s)
    int widest=0;
    {
      std::vector<char> seen((size_t)W*H,0);
      const auto isflat=[&](int x,int y){ if(is_ocean(x,y))return false; float e=dem(x,y);
        for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++){if(!dx&&!dy)continue;int nx=x+dx,ny=y+dy;
          if(dem.inGrid(nx,ny)&&dem(nx,ny)<e)return false;} return true; };
      for(int sy=0;sy<H;sy++)for(int sx=0;sx<W;sx++){
        if(seen[(size_t)sy*W+sx]||!isflat(sx,sy))continue; float e=dem(sx,sy);
        std::vector<std::pair<int,int>> comp; std::vector<std::pair<int,int>> st{{sx,sy}}; seen[(size_t)sy*W+sx]=1;
        int mnx=sx,mxx=sx; bool straddle=false;
        while(!st.empty()){auto[cx,cy]=st.back();st.pop_back();comp.push_back({cx,cy});mnx=std::min(mnx,cx);mxx=std::max(mxx,cx);
          for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++){if(!dx&&!dy)continue;int nx=cx+dx,ny=cy+dy;
            if(!dem.inGrid(nx,ny)||is_ocean(nx,ny)||dem(nx,ny)!=e||!isflat(nx,ny))continue;
            if(!seen[(size_t)ny*W+nx]){seen[(size_t)ny*W+nx]=1;st.push_back({nx,ny});}}}
        if(mnx<s && mxx>=s){ straddle=true; widest=std::max(widest, mxx-mnx+1); }
        (void)straddle;
      }
    }
    // minimal halo for LEFT tile [0,s) and RIGHT tile [s,W)
    auto minhalo=[&](bool left)->int{
      for(int Hh=1; Hh<=W; Hh*=2){
        long conf = left ? owned_matches(ref, replay_window(0, std::min(W,s+Hh)), 0, s)
                         : owned_matches(ref, replay_window(std::max(0,s-Hh), W), s, W);
        if(conf==0) return Hh;
      }
      return -1;
    };
    const int hl=minhalo(true), hr=minhalo(false);
    // "bounded/local" = the owned partition is reproduced with a halo STRICTLY smaller than the full grid
    // (i.e. without seeing the whole domain). max halo tried is < W once s+H first exceeds W, so >=W means
    // only the full-grid window worked -> not local.
    const bool local = (hl>0 && hl<W && hr>0 && hr<W);
    seams++; if(local) bounded++;
    std::cout<<argv[1]<<" seam="<<s<<": widest_straddling_flat="<<widest
             <<"  min_halo(left)="<<hl<<"  min_halo(right)="<<hr
             <<(local?"  BOUNDED-LOCAL":"  NEEDS-FULL-GRID")<<"\n";
  }
  const bool pass = (seams>0 && bounded==seams);
  std::cout<<argv[1]<<": FEASIBILITY "<<(pass?"PASS":"FAIL")<<" ("<<bounded<<"/"<<seams
           <<" seams reproduce serial's owned partition with a bounded local halo)\n";
  return pass ? 0 : 1;
}
