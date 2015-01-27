// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/test/checkadaptation.cc>
#include <dune/grid/../../doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/test/basicunitcube.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

int main (int argc, char *argv[]) try
{

    Dune::GridFactory<FoamGrid<2, 2> > factory;
    BasicUnitCube<2>::insertVertices(factory);
    BasicUnitCube<2>::insertSimplices(factory);

    std::auto_ptr<FoamGrid<2, 2> > grid2d(factory.createGrid());
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("refined0");
    }

    Dune::gridinfo(*grid2d);

    // check grid adaptation interface
    checkAdaptation( *grid2d );

    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("refined-l");
    }
    Dune::gridinfo(*grid2d);
    checkGeometryInFather(*grid2d);
    gridcheck(*grid2d);
    checkIntersectionIterator(*grid2d);

}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
 catch (Exception e) {

    std::cout << e << std::endl;
    return 1;
 }
