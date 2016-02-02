#install headers
install(FILES foamgrid.hh
              dgffoam.cc
              dgffoam.hh
        DESTINATION include/dune/foamgrid)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(foamgrid)
