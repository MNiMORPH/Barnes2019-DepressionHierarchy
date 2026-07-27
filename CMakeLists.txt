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
