// PROOF / regression guard for the flat-label PARTITION rule (SPLIT_INVARIANT_FLATS_PLAN.md ROAD A).
// A full pit-up replay of the flood's cell-labelling, reproducing serial's radix order EXACTLY, compared
// partition-wise (namespace-free) to serial over EVERY land cell. If this ever prints PARTITION-DIFFER the
// flood/radix order changed and dephier_stitch's DH_FLAT_PARTITION_REPLAY pass would no longer reproduce
// serial. Proven on 49 DEMs (24 kerry/testdem fixtures + Corsica + adversarial fractals). Exits nonzero on
// any PARTITION-DIFFER so it can gate as a CTest.
//
// Faithful to dephier.hpp FloodAndAssignDepressions (lines 360-536) + radix_heap.hpp (391-399):
//  * seeds pushed at init: every OCEAN cell, and every land cell with NO strictly-lower neighbour (a
//    land_seed / pit -- this INCLUDES all is_flat interiors), each at its dem elevation.
//  * process buckets by elevation ASCENDING; within a bucket sort ASC by cell index and consume pop_back
//    (HIGHEST-INDEX first); a cell labelled into the CURRENT elevation is appended and popped LIFO.
//  * a popped NO_DEP cell becomes a NEW depression (a pit); a popped labelled/ocean cell propagates its
//    label to NO_DEP neighbours. (Outlet discovery is irrelevant to the label partition, so it's dropped.)
// If the replay's partition equals serial's over all land cells, the flat-partition rule + dynamics are
// fully understood -- no confound.
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
  if(argc!=3){ std::cerr<<"usage: "<<argv[0]<<" <dem> <ocean_level>\n"; return 2; }
  rd::Array2D<float> dem(argv[1]);
  const float ocean_level = std::stod(argv[2]);
  const int W=dem.width(), H=dem.height();
  const auto is_ocean=[&](int x,int y){ return dem.isNoData(x,y) || dem(x,y)==ocean_level; };

  // Serial ground truth: leaf-label grid from the flood.
  rd::Array2D<dh_label_t> s_label(W,H,NO_DEP);
  for(int y=0;y<H;y++) for(int x=0;x<W;x++) if(is_ocean(x,y)) s_label(x,y)=OCEAN;
  rd::Array2D<int8_t> s_fd(W,H,rd::NO_FLOW);
  auto S = dh::GetDepressionHierarchy<float,rd::Topology::D8>(dem, s_label, s_fd);

  // ---- full pit-up replay ----
  rd::Array2D<dh_label_t> r(W,H,NO_DEP);
  std::map<float,std::vector<long>> buckets;   // elevation -> cells (min-key ordered)
  const auto push=[&](int x,int y){ buckets[dem(x,y)].push_back(dem.xyToI(x,y)); };
  for(int y=0;y<H;y++) for(int x=0;x<W;x++){
    if(is_ocean(x,y)){ r(x,y)=OCEAN; push(x,y); continue; }
    const float e=dem(x,y); bool has_lower=false;
    for(int dy=-1;dy<=1&&!has_lower;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=x+dx,ny=y+dy;
      if(dem.inGrid(nx,ny)&&dem(nx,ny)<e){ has_lower=true; break; } }
    if(!has_lower) push(x,y);                  // land_seed (pit / flat interior)
  }
  dh_label_t next_label = 1;                    // ocean is 0; fresh basins from 1 (namespace-free compare)
  while(!buckets.empty()){
    auto it=buckets.begin(); const float e=it->first;
    std::vector<long> cur=std::move(it->second); buckets.erase(it);
    std::sort(cur.begin(),cur.end());          // ascending index; pop_back = highest first
    while(!cur.empty()){
      const long ci=cur.back(); cur.pop_back();
      int cx,cy; dem.iToxy(ci,cx,cy);
      dh_label_t cl=r(cx,cy);
      if(cl==NO_DEP){ cl=next_label++; r(cx,cy)=cl; }   // popped NO_DEP -> new pit
      for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=cx+dx,ny=cy+dy;
        if(!dem.inGrid(nx,ny) || r(nx,ny)!=NO_DEP) continue;
        r(nx,ny)=cl; const float ne=dem(nx,ny);
        if(ne==e) cur.push_back(dem.xyToI(nx,ny)); else buckets[ne].push_back(dem.xyToI(nx,ny)); }
    }
  }

  // ---- partition compare over ALL land cells (namespace-free bijection s_label <-> r) ----
  std::map<dh_label_t,dh_label_t> s2r, r2s;
  long land=0, merge=0, split=0;
  for(int y=0;y<H;y++) for(int x=0;x<W;x++){
    if(is_ocean(x,y)) continue; land++;
    const dh_label_t s=s_label(x,y), rr=r(x,y);
    // OCEAN vs basin is a real difference; count it in both directions.
    auto a=s2r.find(s); if(a==s2r.end()) s2r[s]=rr; else if(a->second!=rr) merge++;
    auto b=r2s.find(rr); if(b==r2s.end()) r2s[rr]=s; else if(b->second!=s) split++;
  }
  const bool ok = (merge+split==0);
  std::cout<<argv[1]<<": land="<<land<<" partition serial-merges(replay-splits)="<<merge
           <<" replay-merges(serial-splits)="<<split
           <<(ok?"  PARTITION-MATCH":"  PARTITION-DIFFER")<<"\n";
  int shown=0;
  for(int y=0;y<H && shown<8;y++) for(int x=0;x<W && shown<8;x++){
    if(is_ocean(x,y)) continue; const dh_label_t s=s_label(x,y), rr=r(x,y);
    if(s2r[s]!=rr || r2s[rr]!=s){
      std::cout<<"  ("<<x<<","<<y<<") elev="<<dem(x,y)<<" serial="<<(long)s<<" replay="<<(long)rr<<"\n"; shown++; }
  }
  return ok ? 0 : 1;
}
