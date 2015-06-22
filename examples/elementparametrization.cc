/** construct a grid with vertices on the unit circle and element parametrization */

#include <config.h>
#include <iostream>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/vtk.hh>

#include <dune/foamgrid/foamgrid.hh>

/**
 * \brief Mapping class mapping from a secant on the unit circle onto the circle
 */
template<typename ctype, int dimgrid, int dimworld>
class UnitCircleMapping :
  public Dune::VirtualFunction<Dune::FieldVector<ctype, dimgrid>, Dune::FieldVector<ctype, dimworld> >
{
  double fromAngle_;
  double toAngle_;

public:
  UnitCircleMapping(double fromAngle, double toAngle) : fromAngle_(fromAngle), toAngle_(toAngle) {}

  ~UnitCircleMapping() {}
  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  void evaluate(const Dune::FieldVector<ctype,dimgrid>& x, Dune::FieldVector<ctype,dimworld>& y) const
  {
    double angle = fromAngle_ + x[0]*(toAngle_ - fromAngle_);
    y = {std::cos(angle), std::sin(angle)};
  }

};

/**
 * \brief Mapping class mapping from a triangle with points on the unit sphere onto the sphere
  *       with theta in [0, pi] and phi in [0, 2*pi)
 */
template<typename ctype, int dimgrid, int dimworld>
class UnitSphereMapping :
  public Dune::VirtualFunction<Dune::FieldVector<ctype, dimgrid>, Dune::FieldVector<ctype, dimworld> >
{
  const std::array<Dune::FieldVector<double, 3>, 3> vertices_;

public:
  UnitSphereMapping(std::array<Dune::FieldVector<double, 3>, 3> vertices) : vertices_(vertices) {}

  ~UnitSphereMapping() {}
  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  void evaluate(const Dune::FieldVector<ctype,dimgrid>& x, Dune::FieldVector<ctype,dimworld>& y) const
  {
    // calculate global coordinate
    Dune::FieldVector<double, 3> shapeFunctions = evaluateShapeFunctions(x);
    y = {0,0,0};
    for(size_t i = 0; i < y.size(); i++)
      for(size_t j = 0; j < 3; j++)
        y[j] += vertices_[i][j]*shapeFunctions[i];
    // project it on the unit sphere
    y /= y.two_norm();
  }

  inline Dune::FieldVector<double, 3> evaluateShapeFunctions(const Dune::FieldVector<ctype,dimgrid>& x) const
  {
    Dune::FieldVector<double, 3> out;
    out[0] = 1.0;
    for (size_t i=0; i<2; i++)
    {
      out[0]  -= x[i];
      out[i+1] = x[i];
    }
    return out;
  }

};

/**
 * \brief Method to calculate the vector update for a single time step advance
 */
template<class GridView, class Mapper>
void evolve (const GridView& gridView,
           const Mapper& mapper,
           std::vector<double>& temperature,
           const double lambda,
           double& dt)
{
  // allocate a temporary vector for the update
  std::vector<double> update(temperature.size());
  std::fill(update.begin(), update.end(), 0.0);

  // initialize dt very large
  dt = std::numeric_limits<double>::max();
  double h = std::numeric_limits<double>::max();

  for (auto&& element : elements(gridView))
  {
    int eIdx = mapper.index(element);
    // iterator over all intersections
    for (auto&& is : intersections(gridView, element))
    {
      // index of the neighbour
      int nIdx = mapper.index(is.outside());

      // calculate distance between the midpoints
      auto eCenter = element.geometry().center();
      auto nCenter = is.outside().geometry().center();
      auto isCenter = is.geometry().center();
      eCenter -= isCenter; nCenter -= isCenter;
      double dist = eCenter.two_norm() + nCenter.two_norm();

      //approximate h as the distance to the neihgbour center
      h = std::min(h, dist);

      // approximate gradient
      double gradTn = (temperature[nIdx] - temperature[eIdx])/dist;

      // add to update
      update[eIdx] += lambda*gradTn;
    }
  }

  // CFL criterion
  dt = std::min(dt, h*h/2.0/lambda);

  // scale dt with safety factor
  dt *= 0.99;

  // update the concentration vector
  for (unsigned int i=0; i<temperature.size(); ++i)
    temperature[i] += dt*update[i];
}

int main (int argc, char *argv[]) try
{
  typedef Dune::FoamGrid<1, 2> Grid;
  typedef Grid::ctype ctype;
  const int dimgrid = Grid::dimension;
  const int dimworld = Grid::dimensionworld;

  // Start grid creation
  Dune::GridFactory<Grid> factory;

  // The list of grid vertex positions
  int numVertices = 3;
  std::vector<Dune::FieldVector<double, dimworld> > vertices({{0.0, 1.0},
                                                              {-0.5*std::sqrt(3), -0.5},
                                                              {0.5*std::sqrt(3), -0.5}});

  // Create the grid vertices
  for (int i=0; i<numVertices; i++)
    factory.insertVertex(vertices[i]);

  // Create the element geometries
  Dune::GeometryType type(1);
  int numElements = 3;
  std::vector<std::vector<unsigned int> > cornerIDs = {{0,1},{1,2},{2,0}};

  double angle = M_PI/2;
  for(int i = 0; i < numElements; i++)
  {
    factory.insertElement(type,
                         cornerIDs[i],
                         std::shared_ptr<Dune::VirtualFunction<Dune::FieldVector<ctype, dimgrid>, Dune::FieldVector<ctype, dimworld> > >
                           (new UnitCircleMapping<ctype, dimgrid, dimworld>(angle, angle+2.0/3.0*M_PI)));
    angle += 2.0/3.0*M_PI;
  }

  // create the grid
  auto grid = factory.createGrid();

  // output VTK
  Dune::VTKWriter<Grid::LeafGridView > writer(grid->leafGridView(), Dune::VTK::nonconforming);
  writer.write("initial");

  // refine the grid
  grid->globalRefine(1);

  // output VTK
  writer.write("refine-1");

  // coarden and then refine four times
  grid->globalRefine(-1);
  grid->globalRefine(4);

  // output VT
  writer.write("refine-4");

  // Solve the heat equation on the refined grid
  // dT/dt + div(lambda*grad(T)) = 0
  // using a finite volume method with an explicit Euler time discretization

  // make a mapper for codim 0 entities in the leaf grid
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGElementLayout>
  mapper(*grid);

  // the primary variable vector
  std::vector<double> temperature(mapper.size());

  // initial conditions
  temperature[0] = 1.0;

  // the time
  double t = 0.0;
  double dt;
  double tEnd = 1.0;
  int timestep = 0;

  // write output only every nth timestep
  int episode = 10;

  // the heat conductivity
  double lambda = 1.0;

  // Write pvd header
  Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView(), "temperature", ".", "");
  vtkWriter.addCellData(temperature, "celldata");
  vtkWriter.write(t);

  // do the time integration
  while(t <= tEnd)
  {
    // apply finite volume scheme
    evolve(grid->leafGridView(), mapper, temperature, lambda, dt);

    //one time step forward
    t += dt;

    //write a vtk
    if(!(timestep%episode))
      vtkWriter.write(t);

    // Output some infos
    std::cout << "Time step " << timestep << " done, t = " << t << ", dt = " << dt << std::endl;
    // Increment time step counter
    ++timestep;
  }

}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {
  std::cout << e << std::endl;
  return 1;
}
