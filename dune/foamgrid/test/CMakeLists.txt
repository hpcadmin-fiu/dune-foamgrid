add_executable(foamgrid-test EXCLUDE_FROM_ALL foamgrid-test.cc)
add_executable(global-refine-test EXCLUDE_FROM_ALL global-refine-test.cc)
add_executable(local-refine-test EXCLUDE_FROM_ALL local-refine-test.cc)

target_link_libraries(foamgrid-test ${DUNE_LIBS})
target_link_libraries(global-refine-test ${DUNE_LIBS})
target_link_libraries(local-refine-test ${DUNE_LIBS})

set(CMAKE_BUILD_TYPE Debug)
