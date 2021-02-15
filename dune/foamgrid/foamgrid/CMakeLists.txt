#install headers
install(FILES foamgrid.impl.hh
              foamgridedge.hh
              foamgridelements.hh
              foamgridentity.hh
              foamgridentitypointer.hh
              foamgridentityseed.hh
              foamgridfactory.hh
              foamgridgeometry.hh
              foamgridhierarchiciterator.hh
              foamgridindexsets.hh
              foamgridintersectioniterators.hh
              foamgridintersections.hh
              foamgridleafiterator.hh
              foamgridleveliterator.hh
              foamgridnulliteratorfactory.hh
              foamgridvertex.hh
         DESTINATION include/dune/foamgrid/foamgrid)


target_sources(dunefoamgrid PRIVATE
     foamgrid.cc
)