// flat_relabel_prototype -- THROWAWAY (ENH-7 prototype). Measures how much a GEOMETRY-DETERMINISTIC
// flat-label rule would change the serial DepressionHierarchy's per-cell labels.
//
// The rule (split-invariant by construction -- pure geometry, no flood order): every FLAT cell (no
// strictly-lower D8 neighbour) is assigned to the depression PIT that is nearest over same-elevation D8
// adjacency (multi-source BFS from all pit cells), breaking ties by the lower pit-cell index. Non-flat
// cells keep their steepest-descent (already deterministic) label. This is the label analog of the
// resolve_flats flow rule: which pit *claims* a flat cell, decided by distance-within-the-flat instead of
// which flood front happened to arrive first.
//
// Reported number = flat cells whose deterministic owner differs from serial's flood label. That is the
// "how much does adopting this rule change serial output?" cost that decides whether ENH-7 is worth it.
// (A tiled build computes the SAME BFS via seam-exchanged relaxations -- like ENH-1 -- so it reproduces
// this deterministic label set exactly; the rule is split-invariant, so distributed == serial under it.)
#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <climits>
#include <iostream>
#include <queue>
#include <tuple>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;
using dh::OCEAN; using dh::NO_DEP;

static rd::Array2D<dh::dh_label_t> ocean_labels(const rd::Array2D<float>& dem, float ol){
  rd::Array2D<dh::dh_label_t> L(dem.width(), dem.height(), NO_DEP);
  for(unsigned i=0;i<dem.size();i++) if(dem.isNoData(i) || dem(i)==ol) L(i)=OCEAN;
  return L;
}

int main(int argc, char** argv){
  if(argc<3){ std::cerr<<"usage: "<<argv[0]<<" <dem> <ocean_level>\n"; return 1; }
  rd::Array2D<float> dem(argv[1]);
  const float ol = std::stof(argv[2]);
  const int W=dem.width(), H=dem.height();

  auto lab = ocean_labels(dem, ol);
  rd::Array2D<int8_t> fd(W,H,rd::NO_FLOW);
  auto G = dh::GetDepressionHierarchy<float,rd::Topology::D8>(dem, lab, fd);  // lab := flood labels

  const auto is_flat=[&](int x,int y)->bool{
    if(dem.isNoData(x,y)) return false;
    const float e=dem(x,y);
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=x+dx,ny=y+dy;
      if(dem.inGrid(nx,ny) && !dem.isNoData(nx,ny) && dem(nx,ny)<e) return false; }
    return true;                                            // no strictly-lower neighbour
  };
  const auto floodpit=[&](int x,int y)->long{
    const auto c=lab(x,y); if(c==OCEAN || c>=G.size()) return -1;
    const auto pc=G[c].pit_cell; return pc==dh::NO_VALUE ? -2 : (long)pc;
  };

  // Multi-source Dijkstra over same-elevation flat adjacency, key = (dist, pit-cell-index).
  const long BIG = (long)W*H + 1;
  std::vector<long> bestkey(W*H, LONG_MAX), bestpit(W*H, -1);
  using Item = std::tuple<int,long,int,int>;               // dist, pit, x, y
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
  for(dh::dh_label_t i=1;i<G.size();i++){
    const auto pc=G[i].pit_cell; if(pc==dh::NO_VALUE) continue;
    int px,py; dem.iToxy(pc,px,py); if(!is_flat(px,py)) continue;
    pq.push({0,(long)pc,px,py});
  }
  const auto key=[&](int d,long pit){ return (long)d*BIG + pit; };
  while(!pq.empty()){
    auto [d,pit,x,y]=pq.top(); pq.pop();
    const int idx=y*W+x; const long k=key(d,pit);
    if(k>=bestkey[idx]) continue;
    bestkey[idx]=k; bestpit[idx]=pit;
    const float e=dem(x,y);
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ if(!dx&&!dy)continue; int nx=x+dx,ny=y+dy;
      if(!dem.inGrid(nx,ny)||dem.isNoData(nx,ny)) continue;
      if(dem(nx,ny)!=e || !is_flat(nx,ny)) continue;
      if(key(d+1,pit) < bestkey[ny*W+nx]) pq.push({d+1,pit,nx,ny});
    }
  }

  long flat=0, covered=0, changed=0, uncovered=0;
  for(int y=0;y<H;y++) for(int x=0;x<W;x++){
    if(!is_flat(x,y) || lab(x,y)==OCEAN) continue;
    flat++;
    const long det=bestpit[y*W+x], fl=floodpit(x,y);
    if(det==-1){ uncovered++; continue; }                  // flat cell no pit-BFS reached (e.g. a sill)
    covered++;
    if(det!=fl) changed++;
  }
  std::cout<<argv[1]<<"  flat_cells="<<flat<<"  bfs_covered="<<covered<<"  uncovered(sills)="<<uncovered
           <<"  deterministic_changes_vs_serial="<<changed<<"\n";
  return 0;
}
