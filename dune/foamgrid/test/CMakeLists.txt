add_definitions(-DDUNE_FOAMGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${DUNE_GRID_EXAMPLE_GRIDS_PATH}\")

add_executable(foamgrid-test foamgrid-test.cc)
add_executable(global-refine-test global-refine-test.cc)
add_executable(local-refine-test local-refine-test.cc)
add_executable(growth-test growth-test.cc)

target_link_libraries(foamgrid-test ${DUNE_LIBS})
target_link_libraries(global-refine-test ${DUNE_LIBS})
target_link_libraries(local-refine-test ${DUNE_LIBS})
target_link_libraries(growth-test ${DUNE_LIBS})

add_test(foamgrid-test foamgrid-test)
add_test(global-refine-test global-refine-test)
add_test(local-refine-test local-refine-test)
add_test(growth-test growth-test)
