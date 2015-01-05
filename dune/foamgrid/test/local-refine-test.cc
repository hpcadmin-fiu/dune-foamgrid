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
    const std::string dune_foamgrid_path = "/home/timokoch/dumux/dune-foamgrid/doc/grids/gmsh/";

    Dune::GridFactory<FoamGrid<2, 2> > factory;
    BasicUnitCube<2>::insertVertices(factory);
    BasicUnitCube<2>::insertSimplices(factory);

    std::cout << "Checking FoamGrid<2, 2> (2d in 2d grid)" << std::endl;
    std::cout << "  Creating grid" << std::endl;
    std::auto_ptr<FoamGrid<2, 2> > grid2d(factory.createGrid());
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView > writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined0");
    }
    gridcheck(*grid2d);
    Dune::gridinfo(*grid2d);

    // check grid adaptation interface
    checkAdaptation(*grid2d);
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView > writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined-l");
    }
    Dune::gridinfo(*grid2d);
    checkGeometryInFather(*grid2d);
    gridcheck(*grid2d);
    checkIntersectionIterator(*grid2d);

    std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;
    std::cout << "  Creating grid" << std::endl;
    std::shared_ptr<FoamGrid<1, 3> > grid1d( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d3d.msh", /*verbose*/ true, false ) );
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined0");
    }
    gridcheck(*grid1d);
    Dune::gridinfo(*grid1d);

    // check grid adaptation interface
    checkAdaptation(*grid1d);
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined-l");
    }
    Dune::gridinfo(*grid1d);
    gridcheck(*grid1d);
    checkGeometryInFather(*grid1d);
    checkIntersectionIterator(*grid1d);
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
 catch (Exception e) {

    std::cout << e << std::endl;
    return 1;
 }
