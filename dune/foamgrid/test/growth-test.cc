// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

template <class Grid>
void checkGridGrowth(Grid& grid)
{
	typedef typename Grid::template Codim<0>::LeafIterator ElementIterator;
    typedef typename ElementIterator::Entity EntityType;
    enum { dimworld = Grid::dimensionworld };
    enum { dim = Grid::dimension };

    const ElementIterator eEndIt = grid.leafGridView().template end<0>();
    for (ElementIterator eIt = grid.leafGridView().template begin<0>(); eIt != eEndIt; ++eIt )
    {
      const EntityType& entity = *eIt;
      auto geo = entity.geometry();
      double eps = 0.0;
      if(geo.center()[1] > eps)
      {
      	Dune::FieldVector<double, dimworld> spawnPoint = geo.center();
      	Dune::FieldVector<double, dimworld> u(0.0);
      	for (int dimIdx = 0; dimIdx < dimworld; ++dimIdx)
      		u[dimIdx] -= 0.2;
      	spawnPoint += u;
      	grid.mark(1, entity, spawnPoint);
      }
	}

	grid.preGrow();
	grid.grow();
	grid.postGrow();
}

using namespace Dune;

int main (int argc, char *argv[]) 
{
  try
  {
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
	std::cout << "Creating a FoamGrid<1, 2> (1d in 2d grid)" << std::endl;
	std::shared_ptr<FoamGrid<1, 2> > grid1d( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );
	{
        Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("before_growth");
    }
    //gridcheck(*grid1d);
    Dune::gridinfo(*grid1d);
    checkGridGrowth(*grid1d);
    {
        Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("after_growth");
    }
    Dune::gridinfo(*grid1d);
  }
  // //////////////////////////////////
  //   Error handler
  // /////////////////////////////////
  catch (Exception e) {
    std::cout << e << std::endl;
    return 1;
  }
  return 0;
};
