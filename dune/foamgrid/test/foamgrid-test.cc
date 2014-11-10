#include <config.h>

#include <iostream>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/../../doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/foamgrid/foamgrid.hh>


int main (int argc, char *argv[]) try
{
    // dimworld == 2
    //FoamGrid<2>* grid2d = make2DHybridTestGrid<FoamGrid<2> >();
    // paths to gmsh test files
    const std::string dune_grid_path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    {
        std::cout << "Checking FoamGrid<2, 2> (2d in 2d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 2> > grid2d( GmshReader<FoamGrid<2, 2> >::read(dune_grid_path + "curved2d.msh", true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid2d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid2d);
    }
    {
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)" << std::endl;
        // dimworld == 3
        std::cout << "  Creating grid" << std::endl;
        FoamGrid<2, 3>* grid3d = make2Din3DHybridTestGrid<FoamGrid<2, 3> >();

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid3d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid3d);
    }
    {
        std::cout << "Checking FoamGrid<1, 2> (1d in 2d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 2> > grid12( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid12);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid12);
    }
    {
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > grid13( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d3d.msh", true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid13);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid13);
    }
    {
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)" << std::endl;

        // dimworld == 3,  and a grid containing a T-Junction
        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 3> > gridTJunction( GmshReader<FoamGrid<2, 3> >::read(dune_foamgrid_path + "tjunction2d.msh", true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*gridTJunction);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*gridTJunction);
    }
    {
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > gridStar( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "bifurcation1d3d.msh", true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*gridStar);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*gridStar);
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Exception e) {
    std::cout << e << std::endl;
    return 1;
}
