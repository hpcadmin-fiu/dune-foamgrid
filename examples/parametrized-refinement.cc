// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
 * \brief Construct a FoamGrid with element parametrizations given by a closed-form function
 */

#include <config.h>
#include <iostream>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/foamgrid/foamgrid.hh>

using namespace Dune;

/**
 * \brief Mapping from R^d to the graph of a given function
 */
template<int dim, int dimworld>
class GraphMapping :
  public Dune::VirtualFunction<Dune::FieldVector<double, dim>, Dune::FieldVector<double, dimworld> >
{
  // Element corners in the global parameter domain
  std::array<FieldVector<double,dim>, dim+1> corners_;

  std::function<FieldVector<double,3>(FieldVector<double,2>)> graph_;

public:
  GraphMapping(std::array<FieldVector<double,dim>, dim+1> corners,
               std::function<FieldVector<double,3>(FieldVector<double,2>)> graph)
  : corners_(corners), graph_(graph)
  {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  void evaluate(const Dune::FieldVector<double,dim>& x, Dune::FieldVector<double,dimworld>& y) const
  {
    // Linear interpolation between the corners
    auto globalX = corners_[0];
    for (size_t i=0; i<x.size(); i++)
      for (int j=0; j<dim; j++)
        globalX[j] += x[i]*(corners_[i+1][j]-corners_[0][j]);
    y = graph_(globalX);
  }

};

int main (int argc, char *argv[]) try
{
  const int dim      = 2;
  const int dimworld = 3;

  // Global parametrization function
  auto parametrization = [](const FieldVector<double,2>& x) -> FieldVector<double,3>
                           {return {x[0], x[1], 0.2*exp(-fabs(x.two_norm()))*cos(4.5*M_PI*x.two_norm())};};


  // Create the grid
  typedef Dune::FoamGrid<dim, dimworld> Grid;

  // Start grid creation
  Dune::GridFactory<Grid> factory;

  // The list of grid vertex positions
  FieldVector<double,2> p0 = {-1,-1};
  FieldVector<double,2> p1 = { 0,-1};
  FieldVector<double,2> p2 = { 1,-1};
  FieldVector<double,2> p3 = {-1, 0};
  FieldVector<double,2> p4 = { 0, 0};
  FieldVector<double,2> p5 = { 1, 0};
  FieldVector<double,2> p6 = {-1, 1};
  FieldVector<double,2> p7 = { 0, 1};
  FieldVector<double,2> p8 = { 1, 1};

  factory.insertVertex(parametrization(p0));
  factory.insertVertex(parametrization(p1));
  factory.insertVertex(parametrization(p2));
  factory.insertVertex(parametrization(p3));
  factory.insertVertex(parametrization(p4));
  factory.insertVertex(parametrization(p5));
  factory.insertVertex(parametrization(p6));
  factory.insertVertex(parametrization(p7));
  factory.insertVertex(parametrization(p8));

  // Create the element geometries
  Dune::GeometryType triangle;
  triangle.makeTriangle();

  std::array<FieldVector<double,2>, 3> corners0 = {p0, p1, p4};
  std::array<FieldVector<double,2>, 3> corners1 = {p0, p4, p3};
  std::array<FieldVector<double,2>, 3> corners2 = {p1, p2, p5};
  std::array<FieldVector<double,2>, 3> corners3 = {p1, p5, p4};
  std::array<FieldVector<double,2>, 3> corners4 = {p3, p4, p7};
  std::array<FieldVector<double,2>, 3> corners5 = {p3, p7, p6};
  std::array<FieldVector<double,2>, 3> corners6 = {p4, p5, p8};
  std::array<FieldVector<double,2>, 3> corners7 = {p4, p8, p7};

  factory.insertElement(triangle, {0,1,4}, std::make_shared<GraphMapping<dim, dimworld> >(corners0, parametrization));
  factory.insertElement(triangle, {0,4,3}, std::make_shared<GraphMapping<dim, dimworld> >(corners1, parametrization));
  factory.insertElement(triangle, {1,2,5}, std::make_shared<GraphMapping<dim, dimworld> >(corners2, parametrization));
  factory.insertElement(triangle, {1,5,4}, std::make_shared<GraphMapping<dim, dimworld> >(corners3, parametrization));
  factory.insertElement(triangle, {3,4,7}, std::make_shared<GraphMapping<dim, dimworld> >(corners4, parametrization));
  factory.insertElement(triangle, {3,7,6}, std::make_shared<GraphMapping<dim, dimworld> >(corners5, parametrization));
  factory.insertElement(triangle, {4,5,8}, std::make_shared<GraphMapping<dim, dimworld> >(corners6, parametrization));
  factory.insertElement(triangle, {4,8,7}, std::make_shared<GraphMapping<dim, dimworld> >(corners7, parametrization));

  // create the grid
  auto grid = factory.createGrid();

  // output VTK
  Dune::VTKWriter<Grid::LeafGridView > writer(grid->leafGridView());

  for (int i=0; i<5; i++)
  {
    writer.write("refine-" + std::to_string(i));
    grid->globalRefine(1);
  }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {
  std::cout << e << std::endl;
  return 1;
}
