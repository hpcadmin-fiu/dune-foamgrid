# install headers
file(GLOB DUNE_FOAMGRID_HEADERS *.hh)
install(FILES ${DUNE_FOAMGRID_HEADERS}
        DESTINATION include/dune/foamgrid/foamgrid)

dune_library_add_sources(dunefoamgrid SOURCES foamgrid.cc)
