// dephier_flat_extent -- connected-flat size/extent statistics for a DEM.
//
// The distributed DepressionHierarchy resolves flat flowdirs per-tile with an adaptive
// boundary halo (PARALLEL_DEPHIER_ENGINEERING.md section 6, move 2). That halo grows to
// a flat's cross-seam EXTENT, so the largest flat's bounding-box dimension is the
// worst-case per-rank halo. This tool measures the flat-extent distribution on real
// tiles so the tile size (and whether the strict-O(boundary) iterative-exchange fallback
// is needed) can be chosen from data, the way section 8 sized centralized Phase C.
//
// Uses richdem's Barnes-2014 flat labelling (resolve_flats_barnes) to group connected
// flats, then reports each group's cell count and bounding box. Emits a machine-readable
// FLATSTATS line plus a human summary.

#include <richdem/common/Array2D.hpp>
#include <richdem/flowmet/d8_flowdirs.hpp>
#include <richdem/flats/flat_resolution.hpp>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>

namespace rd = richdem;

int main(int argc, char **argv){
  if(argc!=2){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM>\n";
    std::cout<<"Reports connected-flat size/extent stats (PARALLEL_DEPHIER_ENGINEERING.md section 6).\n";
    return -1;
  }
  rd::Array2D<float>   dem(argv[1]);
  rd::Array2D<int8_t>  fd(dem.width(), dem.height(), rd::NO_FLOW);
  rd::d8_flow_directions(dem, fd);
  rd::Array2D<int32_t> flat_mask, labels;
  rd::resolve_flats_barnes(dem, fd, flat_mask, labels);

  // Per flat label: cell count + bounding box (label 0 / noData = not a flat).
  std::map<int32_t,long> cnt;
  std::map<int32_t,int>  xmin,xmax,ymin,ymax;
  for(int y=0;y<dem.height();y++) for(int x=0;x<dem.width();x++){
    const int32_t L = labels(x,y);
    if(L==labels.noData() || L==0) continue;
    if(!cnt.count(L)){ xmin[L]=xmax[L]=x; ymin[L]=ymax[L]=y; }
    cnt[L]++;
    xmin[L]=std::min(xmin[L],x); xmax[L]=std::max(xmax[L],x);
    ymin[L]=std::min(ymin[L],y); ymax[L]=std::max(ymax[L],y);
  }

  long nflats=cnt.size(), flatcells=0, maxcells=0, maxdim=0;
  int  hist[6]={0,0,0,0,0,0};                     // bbox max-dim buckets: <=4,8,16,32,64,>64
  for(const auto &kv : cnt){
    const int L=kv.first;
    flatcells += kv.second;
    maxcells   = std::max(maxcells, kv.second);
    const int w=xmax[L]-xmin[L]+1, h=ymax[L]-ymin[L]+1, d=std::max(w,h);
    maxdim = std::max(maxdim, (long)d);
    hist[ d<=4?0 : d<=8?1 : d<=16?2 : d<=32?3 : d<=64?4 : 5 ]++;
  }
  const long ncells = (long)dem.width()*dem.height();

  std::cout<<"# FLATSTATS columns: ncells flats flatcells maxflat_cells maxflat_dim "
             "dim_le4 dim_le8 dim_le16 dim_le32 dim_le64 dim_gt64\n";
  std::cout<<"FLATSTATS "<<ncells<<" "<<nflats<<" "<<flatcells<<" "<<maxcells<<" "<<maxdim<<" "
           <<hist[0]<<" "<<hist[1]<<" "<<hist[2]<<" "<<hist[3]<<" "<<hist[4]<<" "<<hist[5]<<"\n";
  std::cerr<<"flats="<<nflats<<"  flat cells="<<flatcells<<" ("<<(100.0*flatcells/ncells)<<"%)"
           <<"  largest flat="<<maxcells<<" cells, "<<maxdim<<" cells across (worst-case halo)\n";
  return 0;
}
