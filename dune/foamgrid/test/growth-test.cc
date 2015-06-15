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
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  for (const auto& element : elements(grid.leafGridView()))
  {
    const auto geo = element.geometry();
    if(geo.center()[1] >= 0.0 && geo.center()[0] < -1.0)
    {
      Dune::FieldVector<double, dimworld> growPoint = geo.center();
      Dune::FieldVector<double, dimworld> u(0.0);
      for (int dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        u[dimIdx] += 2.0;
      growPoint += u;
      grid.markForGrowth(element, 0, growPoint);
    }
	}

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl<< "numBoundarySegments before growth: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

	bool elementsWillVanish = grid.preGrow();
  if(elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  // check mightVanish
  for (const auto& element : elements(grid.levelGridView(0)))
    checkHierarchy(element);

	bool newElementGenerated = grid.grow();
  if(!newElementGenerated)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  grid.postGrow();

  // Loop over all levels except the lowest one
  for (int level = 0 ; level <= grid.maxLevel(); ++level )
  {
      for (const auto& element : elements(grid.levelGridView(level)))
        if(element.isNew())
        DUNE_THROW(InvalidStateException,"After postGrow() was called no entity is new, i.e., isNew() == false");
  }

  std::cout << "Boundary intersections after growth: " << std::endl;
  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after growth: " << numBoundarySegments << std::endl<< std::endl;
}

template <class Grid>
void checkGridElementMerge(Grid& grid)
{
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  // connect all facets with their closest neighbor facet (if there is no connection already which gets checked by the grid)
  for (const auto& element : elements(grid.leafGridView()))
  {
    auto otherElement = element;
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      Dune::FieldVector<double, dimworld> center = intersection.geometry().center();
      double dist = std::numeric_limits<double>::max();
      int facetIndex = intersection.indexInInside();
      int facetIndex2 = 0;
      for (const auto& element2 : elements(grid.leafGridView()))
      {
        if(element2 == element)
          continue;
        // find the closest of all facets
        for (const auto& intersection2 : intersections(grid.leafGridView(), element2))
        {
          Dune::FieldVector<double, dimworld> center2 = intersection2.geometry().center();
          Dune::FieldVector<double, dimworld> diff = center;
          diff -= center2;
          if(diff.two_norm() < dist)
          {
            dist = diff.two_norm();
            facetIndex2 = intersection2.indexInInside();
            otherElement = element2;
          }
        }
      }
      grid.markForMerging(element, facetIndex, otherElement, facetIndex2);
    }
  }

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl<< "numBoundarySegments before merge: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  grid.preGrow();
  grid.grow();
  grid.postGrow();

  std::cout << "Boundary intersections after merge: " << std::endl;

  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after merge: " << numBoundarySegments << std::endl<< std::endl;
}

template <class Grid>
void checkGridElementRemoval(Grid& grid)
{
  using namespace Dune;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  int counter = 0;
  // connect all facets with their closest neighbor facet (if there is no connection already which gets checked by the grid)
  const auto eEndIt = grid.leafGridView().template end<0>();
  for (auto eIt = grid.leafGridView().template begin<0>(); eIt != eEndIt; ++eIt, ++counter)
  {
    // delete the 9th and the 10th element
    if(counter == 9 || counter == 10)
      grid.markForRemoval(*eIt);
  }

  bool elementsWillVanish = grid.preGrow();
  if(!elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl << "numBoundarySegments before removal: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  grid.grow();
  grid.postGrow();

  std::cout << "Boundary intersections after removal: " << std::endl;

  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after removal: " << numBoundarySegments << std::endl << std::endl;
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
}
