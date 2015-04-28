// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include <memory>

#include <dune/common/parallel/mpihelper.hh> // include mpi helper class

#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/foamgrid/foamgrid.hh>

// the initial condition c0
template<int dimworld, class ct>
double c0 (const Dune::FieldVector<ct,dimworld>& x)
{
  assert(dimworld>=3);
  return (x[2]>1.5) ? 1.0 : 0.0;
}

// the boundary condition b on inflow boundary
template<int dimworld, class ct>
double b (const Dune::FieldVector<ct,dimworld>& x, double t)
{
  return 0.0;
}

// the vector field u
template<int dimworld, class ct>
Dune::FieldVector<double,dimworld> u (const Dune::FieldVector<ct,dimworld>& x, double t)
{
  return {0, 0, -1};
}

//! initialize the vector of unknowns with initial value
template<class GridView, class Mapper>
void initialize (const GridView& gridView, const Mapper& mapper, std::vector<double>& c)
{
  // first we extract the dimensions of the grid
  const int dimworld = GridView::dimensionworld;

  // type used for coordinates in the grid
  typedef typename GridView::ctype ct;

  // iterate through leaf grid an evaluate c0 at cell center
  for (auto it = gridView.template begin<0>(); it != gridView.template end<0>(); ++it)
  {
    // get geometry
    auto geo = it->geometry();

    // get global coordinate of cell center
    Dune::FieldVector<ct,dimworld> global = geo.center();

    // initialize cell concentration
    c[mapper.index(*it)] = c0(global);
  }
}

template<class GridView, class Mapper>
void evolve (const GridView& gridView, const Mapper& mapper, std::vector<double>& c, double t, double& dt)
{
  // first we extract the dimensions of the grid
  const int dimworld = GridView::dimensionworld;

  // type used for coordinates in the grid
  typedef typename GridView::ctype ct;

  // allocate a temporary vector for the update
  std::vector<double> update(c.size());
  std::fill(update.begin(), update.end(), 0.0);

  // initialize dt very large
  dt = std::numeric_limits<double>::max();

  // compute update vector and optimum dt in one grid traversal
  auto endit = gridView.template end<0>();
  for (auto it = gridView.template begin<0>(); it!=endit; ++it)
  {
    // cell geometry
    auto geo = it->geometry();

    // cell volume, assume linear map here
    double volume = geo.volume();

    // cell index
    int indexi = mapper.index(*it);

    // variable to compute sum of positive factors
    double sumfactor = 0.0;

    // run through all intersections /*@$\gamma_{ij}$@*/ with neighbors and boundary
    auto isend = gridView.iend(*it);
    for (auto is = gridView.ibegin(*it); is!=isend; ++is)
    {
      // get geometry type of face
      auto igeo = is->geometry();

      // get normal vector scaled with volume
      Dune::FieldVector<ct,dimworld> integrationOuterNormal
        = is->centerUnitOuterNormal();
      integrationOuterNormal *= igeo.volume();

      // center of face in global coordinates
      Dune::FieldVector<ct,dimworld> faceglobal = igeo.center();

      // evaluate velocity at face center
      Dune::FieldVector<double,dimworld> velocity = u(faceglobal,t);

      // compute factor occuring in flux formula /*@$u \cdot \nu_{ij} \cdot |\gamma_{ij}| / |K_i|$@*/
      double factor = velocity*integrationOuterNormal/volume;

      // for time step calculation
      if (factor>=0)
        sumfactor += factor;

      // handle interior face
      if (is->neighbor())
      {
        // access neighbor
        int indexj = mapper.index(*is->outside());

        // compute flux from one side only
        if (indexi<indexj)
        {
          // compute factor in neighbor
          auto nbgeo = is->outside()->geometry();
          double nbvolume = nbgeo.volume();
          double nbfactor = velocity*integrationOuterNormal/nbvolume;

          if (factor<0)                         // inflow
          {
            update[indexi] -= c[indexj]*factor;
            update[indexj] += c[indexj]*nbfactor;
          }
          else                         // outflow
          {
            update[indexi] -= c[indexi]*factor;
            update[indexj] += c[indexi]*nbfactor;
          }
        }
      }

      // handle boundary face
      if (is->boundary())
      {
        if (factor<0)                 // inflow, apply boundary condition
          update[indexi] -= b(faceglobal,t)*factor;
        else                 // outflow
          update[indexi] -= c[indexi]*factor;
      }
    }             // end all intersections

    // compute dt restriction
    dt = std::min(dt,1.0/sumfactor);

  }       // end grid traversal

  // scale dt with safety factor
  dt *= 0.99;

  // update the concentration vector
  for (unsigned int i=0; i<c.size(); ++i)
    c[i] += dt*update[i];
}

//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main (int argc , char ** argv) try
{
  typedef Dune::FoamGrid<2,3> Grid;

  std::shared_ptr<Grid> grid;

  std::string path = "../../examples/data/";
  std::string gridFile = "y-grid.msh";
  grid = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read(path + "/" + gridFile));

  for (int i=0; i<4; i++)
    grid->globalRefine(1);

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGElementLayout>
  mapper(*grid);

  // allocate a vector for the concentration
  std::vector<double> c(mapper.size());

  // initialize concentration with initial values
  initialize(grid->leafGridView(),mapper,c);

  // Write pvd header
  Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView(), "concentration", ".", "");
  vtkWriter.addCellData(c,"celldata");
  vtkWriter.write( 0.0 );

  // now do the time steps
  double t=0,dt;
  double tend = 2.0;
  int k=0;
  const double saveInterval = 0.1;
  double saveStep = 0.1;
  int counter = 1;

  while (t<tend)
  {
    // augment time step counter
    ++k;

    // apply finite volume scheme
    evolve(grid->leafGridView(),mapper,c,t,dt);

    // augment time
    t += dt;

    // check if data should be written
    if (t >= saveStep)
    {
      // write data
      vtkWriter.clear();
      vtkWriter.addCellData(c,"celldata");
      vtkWriter.write( t );

      // increase counter and saveStep for next interval
      saveStep += saveInterval;
      ++counter;
    }

    // print info about time, timestep size and counter
    std::cout << "s=" << grid->size(0)
              << " k=" << k << " t=" << t << " dt=" << dt << std::endl;
  }

  // done
  return 0;
}
catch (std::exception & e) {
  std::cout << "STL ERROR: " << e.what() << std::endl;
  return 1;
}
catch (Dune::Exception & e) {
  std::cout << "DUNE ERROR: " << e.what() << std::endl;
  return 1;
}
catch (...) {
  std::cout << "Unknown ERROR" << std::endl;
  return 1;
}
