// dh_canonical -- ID-independent canonical form and invariants for a
// DepressionHierarchy, plus a label-remap helper.
//
// This is the core of the differential oracle (PARALLEL_DEPHIER_PLAN.md section 6).
// The serial and distributed builds label depressions in different namespaces, so
// the tester cannot compare raw node IDs. Instead we reduce a hierarchy to:
//
//   * canonicalize()  -- an ID-independent structural signature. Two hierarchies
//                        with the same tree shape and per-node attributes produce
//                        the same string, regardless of how nodes are numbered.
//                        Used for exact checks: determinism, 1x1-tiling identity,
//                        and post-collapse bit-identity (section 3.2).
//   * invariants()    -- aggregate quantities (total volume, cell counts, leaf/
//                        meta counts) that must be conserved even where the trees
//                        differ structurally, i.e. the *refinement* case pre-
//                        collapse (section 3.1).
//
// relabel() applies a bijective label permutation, remapping every label-valued
// field. It both (a) lets us test canonicalize()'s ID-independence and (b) is the
// same operation the tile-split needs to move per-tile labels into a global
// namespace (section 4).

#pragma once

#include <dephier/dephier.hpp>

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace dhtest {

using richdem::dephier::Depression;
using richdem::dephier::DepressionHierarchy;
using richdem::dephier::dh_label_t;
using richdem::dephier::NO_VALUE;   // == NO_PARENT

//Quantize a floating elevation/volume to a fixed number of decimals so that
//tiny fp differences between serial and distributed builds compare equal.
//Infinities (ocean pit/out elevations) stringify as "inf"/"-inf".
template<class T>
static std::string quant(T v, int decimals){
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os.precision(decimals);
  os<<v;
  return os.str();
}

//Recursively emit the canonical signature of the subtree rooted at `node`.
//Children (binary lchild/rchild plus ocean_linked) are canonicalized and sorted,
//so sibling order and numbering never affect the result. `visited` counts nodes
//reached, to assert full, acyclic coverage.
template<class T>
static std::string canon_node(
  const DepressionHierarchy<T> &deps,
  const dh_label_t node,
  const int decimals,
  std::vector<char> &seen,
  std::size_t &visited
){
  seen.at(node) = 1;
  visited++;

  const auto &d = deps[node];

  std::vector<std::string> kids;
  if(d.lchild!=NO_VALUE) kids.push_back(canon_node(deps, d.lchild, decimals, seen, visited));
  if(d.rchild!=NO_VALUE) kids.push_back(canon_node(deps, d.rchild, decimals, seen, visited));
  for(const auto c: d.ocean_linked) kids.push_back(canon_node(deps, c, decimals, seen, visited));
  std::sort(kids.begin(), kids.end());

  std::ostringstream os;
  os<<"("<<quant(d.pit_elev,decimals)
    <<","<<quant(d.out_elev,decimals)
    <<","<<d.cell_count
    <<","<<quant(d.dep_vol,decimals);
  for(const auto &k: kids)
    os<<k;
  os<<")";
  return os.str();
}

//ID-independent canonical signature of the whole hierarchy. Roots are nodes that
//are never referenced as a child (the ocean is one). The forest signature is the
//sorted concatenation of the roots' subtree signatures. Throws if traversal does
//not reach every node exactly once (a cycle or a multi-parent node -- a
//malformed hierarchy).
template<class T>
std::string canonicalize(const DepressionHierarchy<T> &deps, int decimals=4){
  const std::size_t N = deps.size();

  //Mark which nodes are referenced as a child of some other node.
  std::vector<char> is_child(N, 0);
  for(std::size_t i=0;i<N;i++){
    const auto &d = deps[i];
    if(d.lchild!=NO_VALUE) is_child.at(d.lchild) = 1;
    if(d.rchild!=NO_VALUE) is_child.at(d.rchild) = 1;
    for(const auto c: d.ocean_linked) is_child.at(c) = 1;
  }

  std::vector<char> seen(N, 0);
  std::size_t visited = 0;
  std::vector<std::string> roots;
  for(std::size_t i=0;i<N;i++)
    if(!is_child[i])
      roots.push_back(canon_node(deps, static_cast<dh_label_t>(i), decimals, seen, visited));
  std::sort(roots.begin(), roots.end());

  if(visited!=N){
    std::ostringstream os;
    os<<"canonicalize: traversal reached "<<visited<<" of "<<N
      <<" nodes (malformed hierarchy: cycle or multi-parent node)";
    throw std::runtime_error(os.str());
  }

  std::ostringstream os;
  for(const auto &r: roots)
    os<<r;
  return os.str();
}

//Aggregate quantities that must be conserved even when tree structure differs
//(the refinement case, section 3.1). Volume/cells are summed over roots so
//double-counting of nested nodes is avoided.
struct Invariants {
  std::size_t n_nodes = 0;
  std::size_t n_leaf  = 0;
  std::size_t n_meta  = 0;
  double      total_dep_vol   = 0;   //summed over root depressions
  uint64_t    total_cell_count = 0;  //summed over root depressions
};

template<class T>
Invariants invariants(const DepressionHierarchy<T> &deps){
  const std::size_t N = deps.size();
  Invariants inv;
  inv.n_nodes = N;

  std::vector<char> is_child(N, 0);
  for(std::size_t i=0;i<N;i++){
    const auto &d = deps[i];
    if(d.lchild!=NO_VALUE) is_child.at(d.lchild) = 1;
    if(d.rchild!=NO_VALUE) is_child.at(d.rchild) = 1;
    for(const auto c: d.ocean_linked) is_child.at(c) = 1;
  }

  for(std::size_t i=1;i<N;i++){   //skip ocean (node 0)
    const auto &d = deps[i];
    if(d.lchild==NO_VALUE && d.rchild==NO_VALUE) inv.n_leaf++;
    else                                         inv.n_meta++;
  }

  //Total volume/cells conserved between serial and distributed builds. A node's
  //dep_vol/cell_count accumulates its BINARY children (lchild/rchild) only --
  //CalculateTotalVolumes does not roll up ocean_linked children, which hold
  //separate volume (a depression spilling from a cliff into another). So the
  //total, counted once each, is the sum over every non-ocean node that is NOT a
  //binary child of another: the ocean-linked nodes plus any non-ocean roots.
  //(The ocean's own dep_vol is 0*inf = nan and is excluded.)
  std::vector<char> is_binchild(N, 0);
  for(std::size_t i=0;i<N;i++){
    const auto &d = deps[i];
    if(d.lchild!=NO_VALUE) is_binchild.at(d.lchild) = 1;
    if(d.rchild!=NO_VALUE) is_binchild.at(d.rchild) = 1;
  }
  for(std::size_t i=1;i<N;i++){
    if(is_binchild[i]) continue;
    inv.total_dep_vol    += deps[i].dep_vol;
    inv.total_cell_count += deps[i].cell_count;
  }

  return inv;
}

//Apply a bijective label permutation perm[old]=new (perm.size()==deps.size()).
//Remaps every label-valued field. Preserves NO_VALUE/NO_PARENT sentinels.
template<class T>
DepressionHierarchy<T> relabel(
  const DepressionHierarchy<T> &deps,
  const std::vector<dh_label_t> &perm
){
  const auto mapv = [&](dh_label_t x){ return x==NO_VALUE ? NO_VALUE : perm[x]; };

  DepressionHierarchy<T> out(deps.size());
  for(std::size_t i=0;i<deps.size();i++){
    Depression<T> d = deps[i];             //copy attributes
    d.parent    = mapv(d.parent);          //NO_PARENT==NO_VALUE, preserved
    d.odep      = mapv(d.odep);
    d.geolink   = mapv(d.geolink);
    d.lchild    = mapv(d.lchild);
    d.rchild    = mapv(d.rchild);
    d.dep_label = perm[d.dep_label];
    for(auto &c: d.ocean_linked)
      c = perm[c];
    out[perm[i]] = d;
  }
  return out;
}

} // namespace dhtest
