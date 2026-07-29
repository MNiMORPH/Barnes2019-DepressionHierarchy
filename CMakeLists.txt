cmake_minimum_required (VERSION 3.9)

project(dephier)

enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(GDAL)

add_subdirectory(submodules/richdem EXCLUDE_FROM_ALL)

include_directories(submodules)

add_library(dephier INTERFACE)
target_include_directories(dephier INTERFACE include/)
target_compile_features(dephier INTERFACE cxx_auto_type cxx_std_17)

if(USE_GDAL AND GDAL_FOUND)
  add_executable(dephier.exe src/main.cpp)
  target_link_libraries(dephier.exe PRIVATE dephier richdem)

  add_executable(dephier_paper_tests.exe src/dephier_paper_tests.cpp)
  target_link_libraries(dephier_paper_tests.exe PRIVATE dephier richdem)
  target_compile_features(dephier_paper_tests.exe PUBLIC cxx_auto_type cxx_std_17)

  # Tooling for the distributed-build effort (see PARALLEL_DEPHIER_PLAN.md).
  add_executable(dephier_stats.exe tools/dephier_stats.cpp)
  target_link_libraries(dephier_stats.exe PRIVATE dephier richdem)

  add_executable(dephier_flat_extent.exe tools/dephier_flat_extent.cpp)
  target_link_libraries(dephier_flat_extent.exe PRIVATE dephier richdem)

  # Proves the label-free three-field reconstruction of Barnes-2014's flat_mask (the equivalence
  # the ENH-1 distributed flat resolution is built on). See the file header + ENHANCEMENTS.md ENH-1.
  add_executable(flat_mask_reconstruct_test.exe tools/flat_mask_reconstruct_test.cpp)
  target_link_libraries(flat_mask_reconstruct_test.exe PRIVATE dephier richdem)
  add_test(NAME flat_mask_reconstruct COMMAND flat_mask_reconstruct_test.exe)

  # ENH-2 pit-only seeding flag (permit_without_baselevel_seed): a seedless bowl-interior tile must
  # throw with the flag off and flood from its pits with it on. Library-level, no GDAL needed.
  add_executable(bowl_interior_seed_test.exe tools/bowl_interior_seed_test.cpp)
  target_link_libraries(bowl_interior_seed_test.exe PRIVATE dephier richdem)
  add_test(NAME bowl_interior_seed COMMAND bowl_interior_seed_test.exe)

  find_package(Threads REQUIRED)
  add_executable(comm_smoke.exe tools/comm_smoke.cpp)
  target_include_directories(comm_smoke.exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(comm_smoke.exe PRIVATE richdem Threads::Threads)
  add_test(NAME comm_thread_smoke COMMAND comm_smoke.exe)

  # Distributed DH on the thread-comm shim, diffed against a full-grid oracle in-binary
  # (PARALLEL_DEPHIER_ENGINEERING.md section 5, build order step 2; local-first design).
  add_executable(dephier_mpi.exe tools/dephier_mpi.cpp)
  target_include_directories(dephier_mpi.exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(dephier_mpi.exe PRIVATE dephier richdem Threads::Threads)

  # Same harness compiled against the real-MPI comm backend (one rank per process; run with
  # `mpirun -np <ntiles>`). Validates the tree on rank 0. Only built when MPI is available.
  find_package(MPI)
  if(MPI_CXX_FOUND)
    add_executable(dephier_mpi_real.exe tools/dephier_mpi.cpp)
    target_include_directories(dephier_mpi_real.exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(dephier_mpi_real.exe PRIVATE DH_USE_MPI)
    target_link_libraries(dephier_mpi_real.exe PRIVATE dephier richdem MPI::MPI_CXX)

    add_executable(comm_mpi_smoke.exe tools/comm_mpi_smoke.cpp)
    target_include_directories(comm_mpi_smoke.exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(comm_mpi_smoke.exe PRIVATE richdem MPI::MPI_CXX)

    # MPI tests launch real processes via mpirun. --oversubscribe lets ranks exceed cores on a
    # dev box; these are local-only (CI has no MPI runtime, like the richdem-dependent tests).
    add_test(NAME comm_mpi_smoke
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 --oversubscribe
                     $<TARGET_FILE:comm_mpi_smoke.exe>)
    add_test(NAME mpi_real_tree_seam_2rank
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 --oversubscribe
                     $<TARGET_FILE:dephier_mpi_real.exe> ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/edge_seam.dem -9999 20)
    add_test(NAME mpi_real_tree_interior_3rank
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 --oversubscribe
                     $<TARGET_FILE:dephier_mpi_real.exe> ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/edge_interior.dem -9999 13,26)
  endif()

  # The distributed pipeline (built incrementally) must be bit-identical to the in-binary
  # oracle at every stage implemented so far: per-rank PhaseAB (BOUNDARY from a 1-column halo)
  # and the namespace remap (global offsets + label grid). dephier_mpi.exe exits 0 iff every
  # implemented stage matches. Same seam/offset/interior fixtures as the stitch, 2 and 3 ranks.
  foreach(mode seam offset interior)
    add_test(NAME mpi_pipeline_${mode}_2rank
             COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/edge_${mode}.dem -9999 20)
    add_test(NAME mpi_pipeline_${mode}_3rank
             COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/edge_${mode}.dem -9999 13,26)
  endforeach()

  add_executable(dephier_selfcheck.exe tools/dephier_selfcheck.cpp)
  target_link_libraries(dephier_selfcheck.exe PRIVATE dephier richdem)

  add_executable(dephier_stitch.exe tools/dephier_stitch.cpp)
  target_link_libraries(dephier_stitch.exe PRIVATE dephier richdem)

  # ROAD A proof: a full pit-up replay of the flood's cell-labelling, compared partition-wise to serial over
  # every land cell (SPLIT_INVARIANT_FLATS_PLAN.md). Exits nonzero on PARTITION-DIFFER (a CTest guard).
  add_executable(flat_partition_replay_proof.exe tools/flat_partition_replay_proof.cpp)
  target_link_libraries(flat_partition_replay_proof.exe PRIVATE dephier richdem)

  # Lever C feasibility guard: a bounded adaptive-halo windowed replay reproduces serial's flat-label
  # partition on a tile's owned cells (the property the distributed flat-label reconciliation relies on).
  add_executable(flat_seam_feasibility.exe tools/flat_seam_feasibility.cpp)
  target_link_libraries(flat_seam_feasibility.exe PRIVATE dephier richdem)

  # Lever C core (in-process de-risk): per-tile adaptive-halo replay + pit-index labelling (no seam stitch)
  # reproduces serial's flat-label partition in the distributed shape; validates the ASSEMBLED partition == serial.
  add_executable(flat_label_distributed.exe tools/flat_label_distributed.cpp)
  target_link_libraries(flat_label_distributed.exe PRIVATE dephier richdem)

  # Boundary-behavior tests (PARALLEL_DEPHIER_PLAN.md section 6, differential oracle).
  # The tiled stitch must reproduce the serial DepressionHierarchy -- tree AND per-cell
  # flowdirs -- bit-identically. The `seam` fixture puts a depression pit straddling the
  # tile seam (the halo-decision stress case); `offset`/`interior` are near-seam and
  # in-tile controls. dephier_stitch.exe exits 0 iff bit-identical.
  foreach(mode seam offset interior)
    add_test(
      NAME    stitch_edge_${mode}
      COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/edge_${mode}.dem -9999 20
    )
  endforeach()

  # Regression fixtures for the seam-artifact collapse + NoData-as-ocean semantics
  # (2026-07-29 findings). These assert the TREE reproduces serial (STITCH-/MPI-*-MATCH)
  # via stdout, NOT the exit code: on these flat-heavy DEMs the arbitrary in-flat flowdir
  # differs (accepted), which the exit code also gates. The MATCH assertion is robust --
  # it says "tiled == serial" without hardcoding node counts, so it survives a DEM regen
  # and bites precisely if the corresponding fix regresses. Each was shown to fail without
  # its fix (revert the fix -> STITCH-DIFFER).
  #
  # testdem8 has NoData sentinel 9, which ocean_labels reads as OCEAN. split 1 exercises the
  # basin<->NoData-ocean outlet (was inf volume before 261bbcd); split 3 the rim-fragment meta
  # dissolve (collapse Pass B2, 2e86ab7).
  add_test(NAME stitch_nodata_ocean_outlet
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/testdem8.dem 0 1)
  set_tests_properties(stitch_nodata_ocean_outlet PROPERTIES PASS_REGULAR_EXPRESSION "STITCH-MATCH ")
  add_test(NAME stitch_rim_fragment
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/testdem8.dem 0 3)
  set_tests_properties(stitch_rim_fragment PROPERTIES PASS_REGULAR_EXPRESSION "STITCH-MATCH ")

  # corsica (GEBCO clip): richdem reads its ocean cells as NoData, so a coastal degenerate pit
  # whose only escape is the NoData-ocean across the seam must still be contracted (is_seam_artifact
  # NoData escape, 4c3edcb). Split-invariant tree at split 78 in both harnesses.
  add_test(NAME stitch_coastal_escape
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/examples/corsica/corsica.tif 0 78)
  set_tests_properties(stitch_coastal_escape PROPERTIES PASS_REGULAR_EXPRESSION "STITCH-MATCH ")
  add_test(NAME mpi_coastal_escape
           COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/examples/corsica/corsica.tif 0 78)
  set_tests_properties(mpi_coastal_escape PROPERTIES PASS_REGULAR_EXPRESSION "MPI-TREE-MATCH ")

  # The node-count decomposition diagnostic must report DECOMP-CORRECT even when the full
  # signature DIFFERs on the accepted tie-break class -- kerry_test4 split 2 is STITCH-DIFFER
  # but node count 6==6. Locks the diagnostic's separation of a real decomposition error from
  # the meta-order tie-break noise.
  add_test(NAME stitch_decomp_diagnostic
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test4.dem 0 2)
  set_tests_properties(stitch_decomp_diagnostic PROPERTIES PASS_REGULAR_EXPRESSION "STITCH-DECOMP-CORRECT ")

  # Split-invariance = the real acceptance metric (serial-free): the tiled build must produce the SAME tree
  # at different tilings. tools/check_split_invariance.sh compares STITCH-SIG across splits. kerry_test7 is
  # split-invariant (it differs from serial only by a deterministic tie-break, yet is identical across
  # tilings); lock that it stays so. The known split-VARIANT cases (kerry_test2/10/11) are the open target.
  add_test(NAME split_invariance_kerry7
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_split_invariance.sh
                   $<TARGET_FILE:dephier_stitch.exe> ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test7.dem 0 2 3 5 7)

  # Volume-correctness bar (VOL-MATCH) for cases that are NOT bit-identical to serial: the tie-break class
  # is STITCH-DIFFER / DECOMP-INCORRECT but MUST stay volume-correct + valid-tree (the acceptance bar).
  # kerry_test11 split 7 is such a case (meta-order tie-break, total_dep_vol 260458 == serial). An invariant
  # lock: it bites if a future change sends the volume to inf (open depression) or a wrong total. The MPI
  # case additionally exercises the distributed Phase D volume reduction (per-rank partials -> rank 0).
  add_test(NAME stitch_tiebreak_volume
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test11.dem 0 7)
  set_tests_properties(stitch_tiebreak_volume PROPERTIES PASS_REGULAR_EXPRESSION "STITCH-VOL-MATCH ")
  add_test(NAME mpi_tiebreak_volume
           COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test11.dem 0 7)
  set_tests_properties(mpi_tiebreak_volume PROPERTIES PASS_REGULAR_EXPRESSION "MPI-VOL-MATCH ")

  # ENH-8: the distributed build reproduces serial's exact tree via the flat-label replay
  # (DH_FLAT_PARTITION_REPLAY: fully per-rank chained-halo pit-index replay + stitch-style canonicalisation,
  # empty-leaf compaction, collapse retired). Three guards, env-gated so the default build is unchanged:
  #  - basic: kerry_test9 split 7 is MPI-TREE-DIFFER without the flag (cell-assignment residual), MATCH with.
  add_test(NAME mpi_flat_partition_replay
           COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test9.dem 0 7)
  set_tests_properties(mpi_flat_partition_replay PROPERTIES
                       ENVIRONMENT "DH_FLAT_PARTITION_REPLAY=1" PASS_REGULAR_EXPRESSION "MPI-TREE-MATCH ")
  #  - chained relay: 5 tiles with a flat wider than a middle tile, so the halo must be relayed THROUGH an
  #    intervening rank (rank 0's halo spans 3 tiles). Bites if the systolic forwarding breaks.
  add_test(NAME mpi_flat_replay_chained
           COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test10.dem 0 2,4,6,8)
  set_tests_properties(mpi_flat_replay_chained PROPERTIES
                       ENVIRONMENT "DH_FLAT_PARTITION_REPLAY=1" PASS_REGULAR_EXPRESSION "MPI-TREE-MATCH ")
  #  - canonicalisation: kerry_test10 split 23 CRASHES under the naive glab(pit) rule (findSet(NO_VALUE));
  #    the stitch-style min-leaf + true-pit stamp makes it MPI-TREE-MATCH. Bites if that rule regresses.
  add_test(NAME mpi_flat_replay_canonical
           COMMAND dephier_mpi.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test10.dem 0 23)
  set_tests_properties(mpi_flat_replay_canonical PROPERTIES
                       ENVIRONMENT "DH_FLAT_PARTITION_REPLAY=1" PASS_REGULAR_EXPRESSION "MPI-TREE-MATCH ")

  # ROAD A proof/guard: the flat-label PARTITION replay must reproduce serial's cell->leaf partition EXACTLY
  # (proven on 49 DEMs). It is the rule dephier_stitch's DH_FLAT_PARTITION_REPLAY relies on; if the flood/
  # radix pop order ever drifts, this bites before the stitch silently stops matching serial. kerry_test10 is
  # a rich multi-basin fixture; the tool exits nonzero on any PARTITION-DIFFER.
  add_test(NAME flat_partition_replay_proof
           COMMAND flat_partition_replay_proof.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test10.dem 0)
  set_tests_properties(flat_partition_replay_proof PROPERTIES PASS_REGULAR_EXPRESSION "PARTITION-MATCH")

  # Lever C feasibility: kerry_test5 has a 16-wide seam-straddling flat; a bounded local halo must reproduce
  # serial's owned partition at every seam (exit 0). Guards that the distributed flat-label reconciliation
  # can stay local (no whole-domain gather).
  add_test(NAME flat_seam_feasibility
           COMMAND flat_seam_feasibility.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test5.dem 0 2 3 5 7)

  # Lever C core: the distributed-shape flat-label algorithm reproduces serial's partition. kerry_test5
  # split 3,6 (3 tiles, 16-wide straddling flat) is a good multi-seam exercise; exit 0 iff DIST-PARTITION-MATCH.
  add_test(NAME flat_label_distributed
           COMMAND flat_label_distributed.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test5.dem 0 3,6)

  # ROAD A integration: DH_FLAT_PARTITION_REPLAY reconstructs serial's exact cell->leaf partition, closing
  # the cell-assignment DIFFER class. kerry_test4 split 5 is STITCH-DIFFER without the flag and STITCH-MATCH
  # with it (bites: disable the replay -> DIFFER). Env-gated so the default build stays byte-identical.
  add_test(NAME stitch_flat_partition_replay
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test4.dem 0 5)
  set_tests_properties(stitch_flat_partition_replay PROPERTIES
                       ENVIRONMENT "DH_FLAT_PARTITION_REPLAY=1"
                       PASS_REGULAR_EXPRESSION "STITCH-MATCH ")

  # Full identity needs BOTH the flat-partition replay AND serial's canonical lowest-out_cell tie-break
  # (dephier.hpp PhaseAB). kerry_test2 split 2 is a meta-vs-ocean_linked tie that stays STITCH-DIFFER with the
  # replay on until serial breaks the out_elev tie by lowest out_cell (shown to bite: with the old
  # flood-first rule it is DIFFER even with the replay). Locks both fixes -- revert either -> STITCH-DIFFER.
  add_test(NAME stitch_replay_outlet_tiebreak
           COMMAND dephier_stitch.exe ${CMAKE_CURRENT_SOURCE_DIR}/test_cases/kerry_test2.dem 0 2)
  set_tests_properties(stitch_replay_outlet_tiebreak PROPERTIES
                       ENVIRONMENT "DH_FLAT_PARTITION_REPLAY=1"
                       PASS_REGULAR_EXPRESSION "STITCH-MATCH ")

  install(TARGETS dephier.exe DESTINATION bin)
  install(TARGETS dephier_paper_tests.exe DESTINATION bin)
  install(TARGETS dephier_stats.exe DESTINATION bin)
else()
  message(WARNING "GDAL is required to build the DepHier example programs; it is not needed to use the DepHier library in another project.")
endif()
