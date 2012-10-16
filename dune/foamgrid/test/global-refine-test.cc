#include <config.h>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/../../doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/test/checkgeometryinfather.cc>

int main (int argc, char *argv[]) try
{
    // dimworld == 2
    //FoamGrid<2>* grid2d = make2DHybridTestGrid<FoamGrid<2> >();
    // path to gmsh test files
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    std::auto_ptr<FoamGrid<2> > grid2d( GmshReader<FoamGrid<2> >::read( path + "curved2d.msh", false, false ) );

    grid2d->globalRefine(1);
    checkGeometryInFather(*grid2d);

    // dimworld == 3
    FoamGrid<3>* grid3d = make2Din3DHybridTestGrid<FoamGrid<3> >();

    grid3d->globalRefine(1);
    checkGeometryInFather(*grid3d);
    
} 
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
 catch (Exception e) {

    std::cout << e << std::endl;
    return 1;
 }
