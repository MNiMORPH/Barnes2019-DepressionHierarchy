// bowl_interior_seed_test -- unit test for the ENH-2 pit-only seeding flag
// (FloodAndAssignDepressions's `permit_without_baselevel_seed`).
//
// A "bowl-interior" tile in a distributed build lies entirely inside a closed basin whose rim is
// beyond the tile on every side: every one of its cells drains INWARD, so it has no OCEAN and no
// BOUNDARY (base-level) seed -- only pit seeds. FloodAndAssignDepressions's exterior-cell check then finds
// exterior_cells==0. Historically that was an unconditional throw ("No OCEAN or BOUNDARY cells
// found"); ENH-2 added `permit_without_baselevel_seed` so such a tile can flood from its pits alone
// (its open top depression is closed later, across the seam). GUARANTEED to occur at 30" scale (any
// endorheic basin larger than a tile swallows one whole), so the distributed build must handle it.
//
// This is a LIBRARY-level test (not the stitch), because the property is a property of the flag: a
// truly seedless label array is hard to manufacture through the tiling harness (a standalone tile's
// grid edge perturbs the seeding), but trivial to hand here. It asserts BOTH directions of the flag,
// so it bites if the flag is removed OR defaulted wrong:
//   * permit=false  -> the exterior_cells==0 throw MUST fire (the pre-ENH-2 behaviour, preserved).
//   * permit=true   -> NO throw; the flood runs from the pit and yields the ocean node + the basin.
#include <dephier/dephier.hpp>
#include <richdem/common/Array2D.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

namespace rd = richdem;
namespace dh = richdem::dephier;

// Build a closed bowl: a high rim around a single low interior pit, with NO ocean labels anywhere
// (every cell NO_DEP). No cell can drain out -> exterior_cells==0 in FloodAndAssignDepressions.
static void make_bowl(rd::Array2D<float> &dem, rd::Array2D<dh::dh_label_t> &label){
  const int W=5, H=5;
  dem   = rd::Array2D<float>(W,H,10.0f);          // rim
  label = rd::Array2D<dh::dh_label_t>(W,H,dh::NO_DEP);  // no OCEAN / BOUNDARY: seedless
  for(int y=1;y<=3;y++) for(int x=1;x<=3;x++) dem(x,y)=3.0f;  // inner shelf
  dem(2,2)=1.0f;                                   // the pit
}

int main(){
  int failures = 0;

  // (1) permit=false -> the seedless tile must throw (pre-ENH-2 behaviour preserved).
  {
    rd::Array2D<float> dem; rd::Array2D<dh::dh_label_t> label;
    make_bowl(dem,label);
    rd::Array2D<int8_t> fd(dem.width(),dem.height(),rd::NO_FLOW);
    dh::DepressionHierarchy<float> deps;
    std::vector<dh::Outlet<float>> outlets;
    bool threw=false;
    try {
      dh::FloodAndAssignDepressions<float,rd::Topology::D8>(dem,label,fd,deps,outlets,/*permit=*/false);
    } catch(const std::runtime_error &){ threw=true; }
    if(!threw){ std::cerr<<"FAIL: permit=false on a seedless bowl did NOT throw\n"; failures++; }
    else       std::cout<<"ok: permit=false throws on a seedless bowl-interior tile\n";
  }

  // (2) permit=true -> no throw; the pit-only flood yields the ocean node (0) + at least one basin.
  {
    rd::Array2D<float> dem; rd::Array2D<dh::dh_label_t> label;
    make_bowl(dem,label);
    rd::Array2D<int8_t> fd(dem.width(),dem.height(),rd::NO_FLOW);
    dh::DepressionHierarchy<float> deps;
    std::vector<dh::Outlet<float>> outlets;
    try {
      dh::FloodAndAssignDepressions<float,rd::Topology::D8>(dem,label,fd,deps,outlets,/*permit=*/true);
    } catch(const std::exception &e){
      std::cerr<<"FAIL: permit=true on a seedless bowl threw: "<<e.what()<<"\n"; failures++;
    }
    // ocean node (0) is always added; the interior pit becomes at least one more depression.
    if(deps.size()<2){ std::cerr<<"FAIL: permit=true produced "<<deps.size()<<" depressions, expected >=2 (ocean + pit)\n"; failures++; }
    else               std::cout<<"ok: permit=true floods from the pit -> "<<deps.size()<<" depressions (ocean + basin)\n";
  }

  if(failures){ std::cerr<<failures<<" failure(s)\n"; return 1; }
  std::cout<<"bowl_interior_seed_test PASSED\n";
  return 0;
}
