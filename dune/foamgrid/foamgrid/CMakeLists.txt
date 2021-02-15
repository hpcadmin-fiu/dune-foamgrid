# install headers
file(GLOB DUNE_FOAMGRID_HEADERS *.hh)
install(FILES ${DUNE_FOAMGRID_HEADERS}
        DESTINATION include/dune/foamgrid/foamgrid)

target_sources(dunefoamgrid PRIVATE foamgrid.cc)
