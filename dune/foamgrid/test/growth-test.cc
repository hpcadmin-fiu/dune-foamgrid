// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include <dune/grid/io/file/gmshreader.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh>

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkadaptation.hh>
#else
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkadaptation.cc>
#endif

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
      growPoint += Dune::FieldVector<double, dimworld>(2.0);

      // compile all vertex indices of the new element
      std::vector<std::size_t> vertices;

      // insert a new vertex
      std::size_t vIdx = grid.insertVertex(growPoint);
      vertices.push_back(vIdx);

      // find the second index
      Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);
      vertices.push_back(mapper.index(element.template subEntity<dim>(0)));

      // insert the new element
      grid.insertElement(Dune::GeometryType(1), vertices);
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
  enum { dim = Grid::dimension };

  // vertex mapper
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);

  // insert an element between the two closest boundary facets
  std::vector<std::size_t> vertices(dim+1);
  double dist = std::numeric_limits<double>::max();
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      auto center = intersection.geometry().center();
      if(!intersection.boundary())
        continue;
      for (const auto& element2 : elements(grid.leafGridView()))
      {
        if(element2 == element)
          continue;

        for (const auto& intersection2 : intersections(grid.leafGridView(), element2))
        {
          if(!intersection2.boundary())
            continue;
          auto center2 = intersection2.geometry().center();
          auto diff = center;
          diff -= center2;
          if(diff.two_norm() < dist)
          {
            dist = diff.two_norm();
            vertices[0] = mapper.index(element.template subEntity<dim>(intersection.indexInInside()));
            vertices[1] = mapper.index(element2.template subEntity<dim>(intersection2.indexInInside()));
          }
        }
      }
    }
  }
  grid.insertElement(Dune::GeometryType(1), vertices);

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl<< "numBoundarySegments before merge: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool newElementInserted = grid.grow();
  if(!newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");
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

  // remove the first inserted element
  auto eIt = grid.leafGridView().template begin<0>();
  grid.removeElement(*eIt);

  bool elementsWillVanish = grid.preGrow();
  if(!elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl << "numBoundarySegments before removal: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool newElementInserted = grid.grow();
  if(newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");
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

template <class Grid>
void checkGridElementGrowthLevel(Grid& grid)
{
  enum { dim = Grid::dimension };

  // first we partially refine the grid at the xMax boundary
  double xMax = std::numeric_limits<double>::min();
  double xMin = std::numeric_limits<double>::max();
  for (const auto& vertex : vertices(grid.leafGridView()))
  {
    xMax = std::max(xMax, vertex.geometry().center()[0]);
    xMin = std::min(xMin, vertex.geometry().center()[0]);
  }

  for (const auto& element : elements(grid.leafGridView()))
    for (const auto& intersection : intersections(grid.leafGridView(), element))
      if(intersection.boundary())
        if(intersection.geometry().center()[0] == xMax)
          grid.mark(1, element);

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();

  std::cout << std::endl << "-------------------------------------------" << std::endl;
  std::cout << "gridinfo after adaptation: " << std::endl;
  Dune::gridinfo(grid);
  std::cout << "-------------------------------------------" << std::endl;

  // vertex mapper
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);
  // Then we do a merge with a level 1 and a level 0 vertex
  std::vector<std::size_t> elementVertices(dim+1);
  for (const auto& element : elements(grid.leafGridView()))
    for (const auto& intersection : intersections(grid.leafGridView(), element))
      if(intersection.boundary())
      {
        if(intersection.geometry().center()[0] == xMax)
        {
          auto vertex = element.template subEntity<dim>(intersection.indexInInside());
          assert(vertex.level() == 1);
          elementVertices[0] = mapper.index(vertex);
        }
        else if(intersection.geometry().center()[0] == xMin)
        {
          auto vertex = element.template subEntity<dim>(intersection.indexInInside());
          assert(vertex.level() == 0);
          elementVertices[1] = mapper.index(vertex);
        }
      }

  grid.insertElement(Dune::GeometryType(1), elementVertices);

  bool elementsWillVanish = grid.preGrow();
  if(elementsWillVanish)
    DUNE_THROW(Dune::InvalidStateException,"grid.preGrow() does not return correct information");

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl << "numBoundarySegments before complex merge: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool newElementInserted = grid.grow();
  if(!newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");

  for (const auto& element : elements(grid.leafGridView()))
    if(element.isNew())
      assert(element.level() == 0);

  grid.postGrow();

  std::cout << "Boundary intersections after complex merge: " << std::endl;

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
  std::cout << "numBoundarySegments after complex merge: " << numBoundarySegments << std::endl << std::endl;
}

using namespace Dune;

int main (int argc, char *argv[])
{
  try
  {
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    std::cout << "Creating a FoamGrid<1, 2> (1d in 2d grid)" << std::endl;
    std::shared_ptr<FoamGrid<1, 2> > grid1d( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

    Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
    writer.write("before_growth");

    // check simple grid growth
    Dune::gridinfo(*grid1d);
    checkGridElementGrowth(*grid1d);
    writer.write("after_growth");
    Dune::gridinfo(*grid1d);

    // check a merger, i.e. inserting an element only with existing vertices
    checkGridElementMerge(*grid1d);
    writer.write("after_merge");
    Dune::gridinfo(*grid1d);

    // check removal of a grid element
    checkGridElementRemoval(*grid1d);
    writer.write("after_removal");
    Dune::gridinfo(*grid1d);

    // check growth when vertices are on different levels
    checkGridElementGrowthLevel(*grid1d);
    writer.write("after_second_growth");
    Dune::gridinfo(*grid1d);

    // do a grid check on a refined grid
    grid1d->globalRefine(4);
    gridcheck(*grid1d);
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
