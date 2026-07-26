// dephier_selfcheck -- validate the canonical-form oracle core (dh_canonical.hpp)
// against real serial output, before it is used to judge the distributed build.
//
// The differential oracle rests on one property: the canonical signature must be
// invariant to how nodes are numbered, yet sensitive to structural change. This
// driver runs the serial GetDepressionHierarchy on a DEM and checks:
//
//   1. Coverage      -- canonicalize() reaches every node exactly once.
//   2. ID-invariance -- a bijective relabel leaves the signature unchanged.
//   3. Sensitivity   -- mutating one node's outlet elevation changes it.
//
// Exit code is nonzero if any check fails.

#include "dh_canonical.hpp"

#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc!=3){
    std::cout<<"Syntax: "<<argv[0]<<" <Input DEM> <Ocean Level>\n";
    return -1;
  }
  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);

  rd::Array2D<float> topo(in_name);
  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP);
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);
  for(unsigned int i=0;i<label.size();i++)
    if(topo.isNoData(i) || topo(i)==ocean_level)
      label(i) = dh::OCEAN;

  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(topo, label, flowdirs);

  int failures = 0;
  const auto check = [&](const char *name, bool ok){
    std::cerr<<(ok?"  PASS  ":"  FAIL  ")<<name<<"\n";
    if(!ok) failures++;
  };

  std::cerr<<"\n";

  // 1. Coverage: canonicalize throws if traversal misses a node.
  std::string sig;
  try {
    sig = dhtest::canonicalize(deps);
    check("coverage: canonical traversal reaches all nodes", true);
  } catch (const std::exception &e){
    std::cerr<<"         "<<e.what()<<"\n";
    check("coverage: canonical traversal reaches all nodes", false);
  }

  // 2. ID-invariance: relabel with a nontrivial bijection (fix ocean at 0,
  //    reverse the rest) and confirm the signature is unchanged.
  {
    const auto N = static_cast<dh::dh_label_t>(deps.size());
    std::vector<dh::dh_label_t> perm(N);
    perm[0] = 0;                                   //keep ocean at 0
    for(dh::dh_label_t i=1;i<N;i++) perm[i] = N-i; //1..N-1 reversed (a bijection)
    const auto relabeled = dhtest::relabel(deps, perm);
    check("id-invariance: relabel leaves signature unchanged",
          dhtest::canonicalize(relabeled)==sig);
  }

  // 3. Sensitivity: bump one leaf's outlet elevation; the signature must change.
  {
    auto mutated = deps;
    bool bumped = false;
    for(std::size_t i=1;i<mutated.size();i++){
      if(mutated[i].lchild==dh::NO_VALUE && mutated[i].rchild==dh::NO_VALUE){
        mutated[i].out_elev += 1.0f;
        bumped = true;
        break;
      }
    }
    if(bumped)
      check("sensitivity: mutating a node changes the signature",
            dhtest::canonicalize(mutated)!=sig);
    else
      std::cerr<<"  SKIP  sensitivity: no leaf depression to mutate\n";
  }

  //Deterministic 64-bit FNV-1a fingerprint of the canonical signature, so the
  //serial output can be regression-checked across code changes (e.g. the
  //Phase B/C refactor) independent of std::hash seeding.
  uint64_t fp = 1469598103934665603ull;
  for(const unsigned char ch: sig){ fp ^= ch; fp *= 1099511628211ull; }

  const auto inv = dhtest::invariants(deps);
  std::cerr<<"\n";
  std::cerr<<"nodes="<<inv.n_nodes<<"  leaves="<<inv.n_leaf<<"  meta="<<inv.n_meta
           <<"  total_cells="<<inv.total_cell_count
           <<"  total_dep_vol="<<inv.total_dep_vol<<"\n";
  std::cerr<<"signature length = "<<sig.size()<<" chars\n";
  std::cout<<"SIGFP "<<in_name<<" "<<sig.size()<<" "<<std::hex<<fp<<std::dec<<"\n";

  std::cerr<<"\n"<<(failures==0 ? "ALL CHECKS PASSED" : "CHECKS FAILED")<<"\n";
  return failures==0 ? 0 : 1;
}
