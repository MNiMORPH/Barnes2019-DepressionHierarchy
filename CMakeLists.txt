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

  install(TARGETS dephier.exe DESTINATION bin)
  install(TARGETS dephier_paper_tests.exe DESTINATION bin)
  install(TARGETS dephier_stats.exe DESTINATION bin)
else()
  message(WARNING "GDAL is required to build the DepHier example programs; it is not needed to use the DepHier library in another project.")
endif()
