#include <config.h>

#include <iostream>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/../../doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/foamgrid/foamgrid.hh>

template<class G>
void traversal (G& grid)
{
  // first we extract the dimensions of the grid
  const int dimgrid = G::dimension;

  // type used for coordinates in the grid
  // such a type is exported by every grid implementation
  typedef typename G::ctype ct;

  // Leaf Traversal
  std::cout << "*** Traverse codim 0 leaves" << std::endl;

  // type of the GridView used for traversal
  // every grid exports a LeafGridView and a LevelGridView
  typedef typename G :: LeafGridView LeafGridView;

  // get the instance of the LeafGridView
  LeafGridView leafView = grid.leafGridView();

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LeafGridView::template Codim<0>::Iterator ElementLeafIterator;

  // iterate through all entities of codim 0 at the leaves
  int count = 0;
  for (ElementLeafIterator it = leafView.template begin<0>();
       it!=leafView.template end<0>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << it->geometry().corner(0)
              << std::endl;
    count++;
  }

  std::cout << "there are/is " << count << " leaf element(s)" << std::endl;

  // Leafwise traversal of codim dim
  std::cout << std::endl;
  std::cout << "*** Traverse codim " << dimgrid << " leaves" << std::endl;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LeafGridView :: template Codim<dimgrid>
  :: Iterator VertexLeafIterator;

  // iterate through all entities of codim 0 on the given level
  count = 0;
  for (VertexLeafIterator it = leafView.template begin<dimgrid>();
       it!=leafView.template end<dimgrid>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    std::cout << "visiting " << gt
              << " at " << it->geometry().corner(0)
              << std::endl;
    count++;
  }
  std::cout << "there are/is " << count << " leaf vertices(s)"
            << std::endl;

  // Levelwise traversal of codim 0
  std::cout << std::endl;
  std::cout << "*** Traverse codim 0 level-wise" << std::endl;

  // type of the GridView used for traversal
  // every grid exports a LeafGridView and a LevelGridView
  typedef typename G :: LevelGridView LevelGridView;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LevelGridView :: template Codim<0>
  :: Iterator ElementLevelIterator;

  // iterate through all entities of codim 0 on the given level
  for (int level=0; level<=grid.maxLevel(); level++)
  {
    // get the instance of the LeafGridView
    LevelGridView levelView = grid.levelGridView(level);

    count = 0;
    for (ElementLevelIterator it = levelView.template begin<0>();
         it!=levelView.template end<0>(); ++it)
    {
      Dune::GeometryType gt = it->type();
      std::cout << "visiting " << gt
                << " with first vertex at " << it->geometry().corner(0)
                << std::endl;
      count++;
    }
    std::cout << "there are/is " << count << " element(s) on level "
              << level << std::endl;
    std::cout << std::endl;
  }
  // Iterate over all intersections
  std::cout << std::endl;
  std::cout << "*** Traverse intersections with level iterator" << std::endl;
  LevelGridView levelView = grid.levelGridView(0);

  typedef typename LevelGridView::IntersectionIterator LevelIntersectionIterator;
  typedef typename LevelGridView::template Codim<0>::Iterator ElementLevelIterator;

  for(ElementLevelIterator it = levelView.template begin<0>();
        it!=levelView.template end<0>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << it->geometry().corner(0)
              << " and second vertex at " << it->geometry().corner(1);
    if(dimgrid==2)
        std::cout << " and third vertex at " << it->geometry().corner(2);
    std::cout << std::endl;


    count = 0;
    for (LevelIntersectionIterator is = levelView.ibegin(*it);
            is!=levelView.iend(*it); ++is)
    {
        if(is->neighbor()) {
            std::cout << "found neighbor with first vertex at: "
                      << is->outside()->geometry().corner(0) << " and second vertex at: "
                      << is->outside()->geometry().corner(1);
            if(dimgrid==2)
                std::cout << " and third vertex at " << is->outside()->geometry().corner(2);
            std::cout << std::endl;
            ++count;
        } else if(is->boundary()) {
            std::cout << "    this is a boundary intersection." << std::endl;
        }
    }
    std::cout << "This element knows about " << count << " neighbors." << std::endl << std::endl;
  }

  // Iterate over all intersections
  std::cout << std::endl;
  std::cout << "*** Traverse intersections with leaf iterator" << std::endl;

  typedef typename LeafGridView::IntersectionIterator LeafIntersectionIterator;
  typedef typename LeafGridView::template Codim<0>::Iterator ElementLeafIterator;

  for(ElementLeafIterator it = leafView.template begin<0>();
        it!=leafView.template end<0>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    std::cout << "visiting leaf " << gt
              << " with first vertex at " << it->geometry().corner(0)
              << " and second vertex at " << it->geometry().corner(1);
    if(dimgrid==2)
        std::cout << " and third vertex at " << it->geometry().corner(2);
    std::cout << std::endl;

    count = 0;
    for (LeafIntersectionIterator is = leafView.ibegin(*it);
            is!=leafView.iend(*it); ++is)
    {
        if(is->neighbor()){
            std::cout << "    found neighbor with first vertex at: "
                      << is->outside()->geometry().corner(0)
                      << " and second vertex at: "
                      << is->outside()->geometry().corner(1);
            if(dimgrid==2)
                std::cout << " and third vertex at " << is->outside()->geometry().corner(2);
            std::cout << std::endl;
            ++count;
        } else if(is->boundary()) {
            std::cout << "    this is a boundary intersection." << std::endl;
        }
    }
    std::cout << "This element knows about " << count << " neighbors." << std::endl << std::endl;
  }

  // //Check IndexSet
  // for (unsigned int codim = 0; codim <= LevelGridView::dimension; ++codim)
  // {
  //       // walk over all geometry types in the codimension
  //       typedef typename LevelGridView::IndexSet::Types GTV;
  //       GTV gtv = leafView.indexSet().types(codim);
  //       for (typename GTV::const_iterator it = gtv.begin(); it != gtv.end(); ++it)
  //       {
  //         std::cout << "gtv[0]= " << gtv[0] << std::endl;
  //       }
  // }


}

int main (int argc, char *argv[]) try
{
    // paths to gmsh test files
    const std::string dune_grid_path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    {
        std::cout << "Checking FoamGrid<2, 2> (2d in 2d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 2> > grid2d( GmshReader<FoamGrid<2, 2> >::read(dune_grid_path + "curved2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid2d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid2d);

        //std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        //traversal(*grid2d);
    }
    {
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        FoamGrid<2, 3>* grid3d = make2Din3DHybridTestGrid<FoamGrid<2, 3> >();

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid3d);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid3d);

        //std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        //traversal(*grid3d);
    }
    {
        std::cout << "Checking FoamGrid<1, 2> (1d in 2d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 2> > grid12( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid12);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid12);

        //std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        //traversal(*grid12);
    }
    {
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > grid13( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d3d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid13);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid13);

        //std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        //traversal(*grid13);
    }
    {
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)" << std::endl;

        // dimworld == 3,  and a grid containing a T-Junction
        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<2, 3> > gridTJunction( GmshReader<FoamGrid<2, 3> >::read(dune_foamgrid_path + "tjunction2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*gridTJunction);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*gridTJunction);

        //std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        //traversal(*gridTJunction);
    }
    {
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > gridStar( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "bifurcation1d3d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*gridStar);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*gridStar);

        std::cout << "  Check if has multiple neighbor functionality" << std::endl;
        traversal(*gridStar);
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Exception e) {
    std::cout << e << std::endl;
    return 1;
}
