// dephier_stats -- report DepressionHierarchy size statistics for a DEM.
//
// Purpose: bound the "reconciliation footprint" the distributed build rests on
// (PARALLEL_DEPHIER_PLAN.md section 8). The tree the serial build returns is the
// O(#depressions) object a distributed build must gather/centralize; this tool
// measures it directly on real or synthetic tiles so we can decide whether a
// centralized Phase C fits on one node.
//
// It runs the *unmodified* serial GetDepressionHierarchy and inspects only its
// return value, so it never perturbs the library code that serves as the oracle.
//
// Emits a machine-readable "DHSTATS" line (easy to sweep/aggregate) plus a
// human-readable block.
//
// Note: #outlets (internal + cross-tile) is *not* recoverable from the returned
// hierarchy -- it lives inside the build and, for cross-tile links, only exists
// once tiling is introduced. This tool reports the depression-tree footprint
// (nodes + bytes); outlet-graph counts come later, from the tiling step.

#include <dephier/dephier.hpp>

#include <richdem/common/Array2D.hpp>

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc!=3){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level>\n";
    std::cout<<"Reports DepressionHierarchy size statistics (see PARALLEL_DEPHIER_PLAN.md section 8).\n";
    return -1;
  }

  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);

  try {

  rd::Array2D<float> topo(in_name);

  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP);
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);

  //Label the ocean cells -- precondition for GetDepressionHierarchy.
  uint64_t ocean_cells = 0;
  for(unsigned int i=0;i<label.size();i++){
    if(topo.isNoData(i) || topo(i)==ocean_level){
      label(i) = dh::OCEAN;
      ocean_cells++;
    }
  }

  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(topo, label, flowdirs);

  //Classify nodes. deps[0] is the ocean; leaves have no children, meta-
  //depressions (internal nodes) have two.
  uint64_t leaves = 0, metas = 0, ocean_linked = 0;
  for(unsigned int i=1;i<deps.size();i++){
    const auto &d = deps[i];
    if(d.lchild==dh::NO_VALUE && d.rchild==dh::NO_VALUE)
      leaves++;
    else
      metas++;
    ocean_linked += d.ocean_linked.size();
  }

  const uint64_t ncells   = static_cast<uint64_t>(topo.width())*topo.height();
  const uint64_t node_sz  = sizeof(dh::Depression<float>);
  const uint64_t tree_bytes = static_cast<uint64_t>(deps.size())*node_sz;

  //Machine-readable line for sweeps: aggregate across tiles by summing columns.
  std::cout<<"# DHSTATS columns: width height ncells datacells oceancells ndep nleaf nmeta oceanlinked node_bytes tree_bytes\n";
  std::cout<<"DHSTATS "
           <<topo.width()      <<" "
           <<topo.height()     <<" "
           <<ncells            <<" "
           <<topo.numDataCells()<<" "
           <<ocean_cells       <<" "
           <<deps.size()       <<" "
           <<leaves            <<" "
           <<metas             <<" "
           <<ocean_linked      <<" "
           <<node_sz           <<" "
           <<tree_bytes        <<"\n";

  //Human-readable summary.
  std::cerr<<"\n";
  std::cerr<<"DEM cells (w*h)      = "<<ncells<<" ("<<topo.width()<<" x "<<topo.height()<<")\n";
  std::cerr<<"Data cells           = "<<topo.numDataCells()<<"\n";
  std::cerr<<"Ocean cells          = "<<ocean_cells<<"\n";
  std::cerr<<"Depressions (nodes)  = "<<deps.size()<<"  [1 ocean + "<<leaves<<" leaves + "<<metas<<" meta]\n";
  std::cerr<<"Ocean-linked links   = "<<ocean_linked<<"\n";
  std::cerr<<"Depression node size = "<<node_sz<<" bytes\n";
  std::cerr<<"Tree footprint       = "<<tree_bytes<<" bytes ("
           <<(tree_bytes/1048576.0)<<" MiB)\n";
  std::cerr<<"Tree bytes / DEM cell= "<<(static_cast<double>(tree_bytes)/ncells)<<"\n";

  } catch (const std::exception &e){
    //E.g. GetDepressionHierarchy throws when the chosen ocean level yields no
    //ocean cells. Skip this tile cleanly rather than abort so sweeps continue.
    std::cerr<<"dephier_stats: skipped '"<<in_name<<"': "<<e.what()<<"\n";
    return 2;
  }

  return 0;
}
