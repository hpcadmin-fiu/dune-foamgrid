// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/test/checkadaptation.hh>

template <class Grid>
void checkGridElementGrowth(Grid& grid)
{
  using namespace Dune;
	typedef typename Grid::template Codim<0>::LeafIterator ElementIterator;
  typedef typename ElementIterator::Entity EntityType;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  const ElementIterator eEndIt = grid.leafGridView().template end<0>();
  for (ElementIterator eIt = grid.leafGridView().template begin<0>(); eIt != eEndIt; ++eIt )
  {
    const EntityType& entity = *eIt;
    auto geo = entity.geometry();
    if(geo.center()[1] >= 0.0 && geo.center()[0] < -1.0)
    {
      Dune::FieldVector<double, dimworld> growPoint = geo.center();
      Dune::FieldVector<double, dimworld> u(0.0);
      for (int dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        u[dimIdx] += 2.0;
      growPoint += u;
      grid.markForGrowth(entity, 0, growPoint);
    }
	}

	bool elementsWillVanish = grid.preGrow();
  if(elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  // check mightVanish
  typedef typename Grid::template Codim<0>::LevelIterator LevelElementIterator;
  for (LevelElementIterator it = grid.levelGridView(0).template begin<0>(); it != grid.levelGridView(0).template end<0>(); ++ it)
    checkHierarchy(*it);

	bool newElementGenerated = grid.grow();
  if(!newElementGenerated)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  grid.postGrow();

  // Loop over all levels except the lowest one
  for (int level = 0 ; level <= grid.maxLevel(); ++level )
  {
    typedef typename Grid::template Codim<0>::LevelIterator LevelElementIterator;
    for (LevelElementIterator it = grid.levelGridView(level).template begin<0>(); it != grid.levelGridView(level).template end<0>(); ++ it)
      if(it->isNew())
        DUNE_THROW(InvalidStateException,"After postGrow() was called no entity is new, i.e., isNew() == false");
  }
}

template <class Grid>
void checkGridElementMerge(Grid& grid)
{
  typedef typename Grid::template Codim<0>::LeafIterator ElementIterator;
  typedef typename ElementIterator::Entity EntityType;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  // connect all facets with their closest neighbor facet (if there is no connection already which gets checked by the grid)
  const ElementIterator eEndIt = grid.leafGridView().template end<0>();
  for (ElementIterator eIt = grid.leafGridView().template begin<0>(); eIt != eEndIt; ++eIt )
  {
    const EntityType& entity = *eIt;
    ElementIterator entity2It = eIt;
    for(auto isIt = grid.leafGridView().ibegin(*eIt); isIt != grid.leafGridView().iend(*eIt); ++isIt)
    {
      Dune::FieldVector<double, dimworld> center = isIt->geometry().center();
      double dist = std::numeric_limits<double>::max();
      int facetIndex = isIt->indexInInside();
      int facetIndex2 = 0;
      for (ElementIterator eIt2 = grid.leafGridView().template begin<0>(); eIt2 != eEndIt; ++eIt2 )
      {
        if(eIt2 == eIt)
          continue;
        // find the closest of all facets
        for(auto isIt2 = grid.leafGridView().ibegin(*eIt2); isIt2 != grid.leafGridView().iend(*eIt2); ++isIt2)
        {
          Dune::FieldVector<double, dimworld> center2 = isIt2->geometry().center();
          Dune::FieldVector<double, dimworld> diff = center;
          diff -= center2;
          if(diff.two_norm() < dist)
          {
            dist = diff.two_norm();
            facetIndex2 = isIt2->indexInInside();
            entity2It = eIt2;
          }
        }
      }
      grid.markForMerging(entity, facetIndex, *entity2It, facetIndex2);
    }
  }

  grid.preGrow();
  grid.grow();
  grid.postGrow();
}

template <class Grid>
void checkGridElementRemoval(Grid& grid)
{
  using namespace Dune;
  typedef typename Grid::template Codim<0>::LeafIterator ElementIterator;
  typedef typename ElementIterator::Entity EntityType;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  int counter = 0;
  // connect all facets with their closest neighbor facet (if there is no connection already which gets checked by the grid)
  const ElementIterator eEndIt = grid.leafGridView().template end<0>();
  for (ElementIterator eIt = grid.leafGridView().template begin<0>(); eIt != eEndIt; ++eIt, ++counter)
  {
    // delete the 9th and the 10th element
    if(counter == 9 || counter == 10)
      grid.markForRemoval(*eIt);
  }

  bool elementsWillVanish = grid.preGrow();
  if(!elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

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
    checkGridElementGrowth(*grid1d);
    {
        Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("after_growth");
    }
    Dune::gridinfo(*grid1d);
    checkGridElementMerge(*grid1d);
    {
        Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("after_merge");
    }
    Dune::gridinfo(*grid1d);
    checkGridElementRemoval(*grid1d);
    {
        Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("after_removal");
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
