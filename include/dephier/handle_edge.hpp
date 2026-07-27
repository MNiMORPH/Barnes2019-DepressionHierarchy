// Cross-tile edge matching for the distributed DepressionHierarchy stitch.
//
// HandleEdge/HandleCorner are lifted from Richard Barnes' parallel priority-flood and
// generalised to a callback:
//   submodules/richdem/programs/parallel_priority_flood/main.cpp  (HandleEdge, HandleCorner)
//
// They pair the 1-cell perimeter strips of two adjacent tiles under the D8 rule -- edge
// cell i in strip A against neighbours {i-1, i, i+1} in strip B -- so that only the strips
// (O(boundary)) cross the seam, never a tile interior. The per-pair "outlet is the higher
// of the two cells, keep the lowest such per depression pair" rule (the serial PhaseB rule)
// is left to a caller-supplied `emit(label_a, elev_a, label_b, elev_b)`, so the same
// primitive feeds either a spill-elevation graph (filling) or an Outlet database (the
// hierarchy). Global label namespacing (label_offset in the original) is applied by the
// caller before/inside emit.

#pragma once

#include <cstddef>
#include <vector>

namespace richdem {
namespace dephier {

// Match two aligned perimeter strips (equal length) across a seam. For each cell i in strip
// A and its D8 neighbours {i-1, i, i+1} in strip B, invoke emit() with the two cells'
// (label, elevation). emit performs the same-depression skip and the min-of-max-per-pair
// bookkeeping (see the callers).
template<class elev_t, class label_t, class Emit>
void HandleEdge(
  const std::vector<elev_t>  &elev_a,
  const std::vector<elev_t>  &elev_b,
  const std::vector<label_t> &label_a,
  const std::vector<label_t> &label_b,
  Emit emit
){
  const std::size_t len = elev_a.size();
  for(std::size_t i=0; i<len; i++)
    for(long ni=(long)i-1; ni<=(long)i+1; ni++){
      if(ni<0 || ni>=(long)len) continue;
      emit(label_a[i], elev_a[i], label_b[(std::size_t)ni], elev_b[(std::size_t)ni]);
    }
}

// The single diagonal cell where four tiles meet (one A cell vs one B cell). Not exercised
// by column-split tiling (it needs 2-D tiles); provided for the 2-D MPI build.
template<class elev_t, class label_t, class Emit>
void HandleCorner(elev_t elev_a, label_t label_a, elev_t elev_b, label_t label_b, Emit emit){
  emit(label_a, elev_a, label_b, elev_b);
}

} // namespace dephier
} // namespace richdem
