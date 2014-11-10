#include <config.h>

#include <iostream>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/../../doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/foamgrid/foamgrid.hh>

template<class G>
void traversal (G& grid)
{
  // first we extract the dimensions of the grid
  const int dim = G::dimension;                        

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
  std::cout << "*** Traverse codim " << dim << " leaves" << std::endl;

  // Get the iterator type
  // Note the use of the typename and template keywords
  typedef typename LeafGridView :: template Codim<dim>
  :: Iterator VertexLeafIterator;                

  // iterate through all entities of codim 0 on the given level
  count = 0;
  for (VertexLeafIterator it = leafView.template begin<dim>(); 
       it!=leafView.template end<dim>(); ++it)
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
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
  typedef typename LevelGridView::template Codim<0>::Iterator ElementLevelIterator;
  
  ElementLevelIterator endlevelit = levelView.template end<0>();
  for(ElementLevelIterator itlevel = levelView.template begin<0>(); itlevel!=endlevelit; ++itlevel)
  {
    Dune::GeometryType gt = itlevel->type();       
    std::cout << "visiting leaf " << gt                
              << " with first vertex at " << itlevel->geometry().corner(0)
              << " and second vertex at " << itlevel->geometry().corner(1)
              << std::endl;
    count = 0;
    LevelIntersectionIterator isend = levelView.iend(*itlevel);
    for (LevelIntersectionIterator is = levelView.ibegin(*itlevel); is!=isend; ++is)
    {
            EntityPointer neighbor = is->outside();
            std::cout << "found neighbor with first vertex at: " 
                      << neighbor->geometry().corner(0) << " and second vertex at: " 
                      << neighbor->geometry().corner(1) << std::endl; 
            ++count;
    }
    std::cout << "This element knows about " << count << " neighbors." << std::endl << std::endl;
  }

  // Iterate over all intersections
  std::cout << std::endl;
  std::cout << "*** Traverse intersections with leaf iterator" << std::endl;
  typedef typename LeafGridView::IntersectionIterator LeafIntersectionIterator;
  typedef typename LeafGridView::template Codim<0>::Iterator ElementLeafIterator;
  
  ElementLeafIterator endleafit = leafView.template end<0>();
  for(ElementLeafIterator itleaf = leafView.template begin<0>(); itleaf!=endleafit; ++itleaf)
  {
    Dune::GeometryType gt = itleaf->type();       
    std::cout << "visiting leaf " << gt                
              << " with first vertex at " << itleaf->geometry().corner(0)
              << " and second vertex at " << itleaf->geometry().corner(1)
              << std::endl;
    count = 0;
    LeafIntersectionIterator isend = leafView.iend(*itleaf);
    for (LeafIntersectionIterator is = leafView.ibegin(*itleaf); is!=isend; ++is)
    {
            EntityPointer neighbor = is->outside();
            std::cout << "found neighbor with first vertex at: " 
                      << neighbor->geometry().corner(0) << " and second vertex at: " 
                      << neighbor->geometry().corner(1) << std::endl; 
            ++count;
    }
    std::cout << "This element knows about " << count << " neighbors." << std::endl << std::endl;
  }

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
    }
    {
        std::cout << "Checking FoamGrid<2, 3> (2d in 3d grid)" << std::endl;

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
        std::shared_ptr<FoamGrid<1, 2> > grid12( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid12);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid12);
    }
    {
        std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;

        std::cout << "  Creating grid" << std::endl;
        std::shared_ptr<FoamGrid<1, 3> > grid13( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d3d.msh", /*verbose*/ true, false ) );

        std::cout << "  Calling gridcheck" << std::endl;
        gridcheck(*grid13);

        std::cout << "  Calling checkIntersectionIterator" << std::endl;
        checkIntersectionIterator(*grid13);
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
