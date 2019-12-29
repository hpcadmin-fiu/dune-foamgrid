#include <config.h>

#include <iostream>
#include <random>
#include <algorithm>

#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/foamgrid/foamgrid.hh>

// Artifical tree growth following the "space colonialization" algorithm of Runions et al (2007)
namespace TreeGrowth {

template<int dimworld>
class Attractor
{
  using Point = Dune::FieldVector<double, dimworld>;
public:
  // initialize the attractor points in a unit cube
  Attractor(std::size_t numAttractors = 1000, double minHeight = 0.1)
  {
    points_.resize(numAttractors);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (auto& att : points_)
      for (int dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        att[dimIdx] = dist(gen);

    // make sure all point are above the minimum height
    for (auto& att : points_)
      if (att[dimworld-1] < minHeight)
        att[dimworld-1] += dist(gen)*(1.0-minHeight) + (minHeight-att[dimworld-1]);

    // project point outside sphere onto sphere
    Point center(0.5);
    for (auto& att : points_)
      if ((att-center).two_norm2() > 0.25)
        att.axpy((att-center).two_norm()-0.5, center-att);
  }

  struct ForceInfo
  {
    Point direction{};
    std::size_t count = 0;
  };

  template<class GridView>
  std::vector<ForceInfo> computeForces(const GridView& gridView, double maxDist)
  {
    std::vector<ForceInfo> attractorForces(gridView.size(1));
    for (const auto& point : points_)
    {
      std::size_t closestVertex = 0;
      double closestVertexDist = std::numeric_limits<double>::max();
      Point direction{};
      for (const auto& vertex : vertices(gridView))
      {
        const auto pos = vertex.geometry().corner(0);
        const auto orientation = (point-pos);
        const auto dist2 = orientation.two_norm2();
        if (dist2 > maxDist*maxDist)
          continue;

        if (dist2 < closestVertexDist)
        {
          closestVertexDist = dist2;
          direction = orientation;
          closestVertex = gridView.indexSet().index(vertex);
        }
      }

      // if we found a closest vertex that is not the root
      if (closestVertex != 0)
      {
        auto& info = attractorForces[closestVertex];
        info.direction.axpy(1.0/direction.two_norm(), direction);
        info.count++;
      }
    }

    // normalize attractor force
    for (const auto& vertex : vertices(gridView))
    {
      auto& info = attractorForces[gridView.indexSet().index(vertex)];
      if (info.count > 0)
        info.direction /= static_cast<double>(info.count);
    }

    return attractorForces;
  }

  template<class GridView>
  void prune(const GridView& gridView, double minDist)
  {
    points_.erase(std::remove_if(points_.begin(), points_.end(),
      [&](const auto& point){
        for (const auto& vertex : vertices(gridView))
        {
          const auto pos = vertex.geometry().corner(0);
          if ((pos-point).two_norm2() < minDist*minDist)
            return true;
        }
        return false;
      }),
      points_.end());
  }

  // check if the attractor is empty
  bool empty() const
  { return points_.empty(); }

  // print for debugging or visualization
  void dumpToFile(const std::string& fileName)
  {
    std::ofstream outfile(fileName, std::ios::out);
    outfile << "x y z\n";
    std::ostream_iterator<Point> it(outfile, "\n");
    std::copy(points_.begin(), points_.end(), it);
  }

private:
  std::vector<Point> points_;
};

template<int dimworld>
class Tree
{
public:
  using Grid = Dune::FoamGrid<1, dimworld>;
  using Point = Dune::FieldVector<double, dimworld>;
  static constexpr int dimensionworld = dimworld;

  // initialize the tree
  Tree(double dx = 0.01)
  : dx_(dx)
  {
      Dune::GridFactory<Grid> factory;
      {
        Point p0(0.5); p0[dimworld-1] = 0.0;
        Point dir(0.0); dir[dimworld-1] = 1.0;
        Point p1 = p0; p1.axpy(dx, dir);
        factory.insertVertex(p0);
        factory.insertVertex(p1);
      }{
        Point p0(0.3); p0[dimworld-1] = 0.0;
        Point dir(0.0); dir[dimworld-1] = 1.0;
        Point p1 = p0; p1.axpy(dx, dir);
        factory.insertVertex(p0);
        factory.insertVertex(p1);
      }
      factory.insertElement(Dune::GeometryTypes::line, std::vector<unsigned int>({0, 1}));
      factory.insertElement(Dune::GeometryTypes::line, std::vector<unsigned int>({2, 3}));
      grid_ = std::unique_ptr<Grid>(factory.createGrid());
  }

  // get a tree view
  typename Grid::LeafGridView gridView() const
  { return grid_->leafGridView(); }

  // the size of the tree
  std::size_t size() const
  { return grid_->leafGridView().size(0); }

  // get access to the tree
  Grid& grid()
  { return *grid_; }

  // get the step size
  double dx() const
  { return dx_; }

private:
  double dx_;
  std::unique_ptr<Grid> grid_;
};

template<class Tree, class Attractor>
void grow(Tree& tree, Attractor& attractor, const Dune::ParameterTree& params)
{
  const auto maxDist = params.get<double>("MaxDist");
  const auto gridView = tree.gridView();
  // get the info on how to grow
  auto forces = attractor.computeForces(gridView, maxDist);

  // grow the tree
  tree.grid().preGrow();
  for (const auto& element : elements(gridView))
  {
    for (int vIdxLocal = 0; vIdxLocal < 2; ++vIdxLocal)
    {
      const auto vertex = element.template subEntity<1>(vIdxLocal);
      const auto vIdxGlobal = gridView.indexSet().index(vertex);
      // grow at this vertex
      if (forces[vIdxGlobal].count > 0)
      {
        auto pos = vertex.geometry().corner(0);
        pos.axpy(tree.dx(), forces[vIdxGlobal].direction);
        const auto newVIdx = tree.grid().insertVertex(pos);
        tree.grid().insertElement(Dune::GeometryTypes::line, std::vector<unsigned int>({vIdxGlobal, newVIdx}));
        forces[vIdxGlobal].count = 0;
      }
    }
  }

  tree.grid().grow();
  tree.grid().postGrow();

  // prune the attractor
  const auto minDist = params.get<double>("MinDist");
  attractor.prune(gridView, minDist);
}

// compute branch radii according to Murray's law with a given exponent (r_parent^n = r_1^n + r_2^n + ...)
template<class GridView>
void computeRadii(const GridView& gridView, std::vector<double>& radius, const std::vector<std::size_t>& roots,
                  double minRadius = 0.001, double exponent = 2.3)
{
  radius.resize(gridView.size(0));
  std::vector<bool> visited(gridView.size(0), false);
  std::vector<typename GridView::template Codim<0>::Entity::EntitySeed> seeds(gridView.size(0));
  for (const auto& element : elements(gridView))
    seeds[gridView.indexSet().index(element)] = element.seed();
  std::stack<typename GridView::template Codim<0>::Entity> elemStack;

  auto computeNeighbors = [&](const auto& element)
  {
    std::array<std::pair<int, Dune::ReservedVector<std::size_t, 8>>, 2> neighbors{};
    for (const auto& intersection : intersections(gridView, element))
    {
      if (intersection.neighbor())
      {
        const auto nIdx = gridView.indexSet().index(intersection.outside());
        neighbors[intersection.indexInInside()].first += !visited[nIdx] ? 1 : 0;
        neighbors[intersection.indexInInside()].second.push_back(nIdx);
      }
    }

    return neighbors;
  };

  auto handleElement = [&](const auto& element)
  {
    const auto neighbors = computeNeighbors(element);
    const auto eIdx = gridView.indexSet().index(element);
    // first compute radius of this element
    const int completedSide = [&](){
      if (neighbors[0].first)
        return 1;
      else if (neighbors[0].second.empty() && std::find(roots.begin(), roots.end(), eIdx) != roots.end())
        return 1;
      else
        return 0;
    }();

    visited[eIdx] = true;
    radius[eIdx] = [&](){
      if (neighbors[completedSide].second.size() == 0)
        return minRadius;
      else if (neighbors[completedSide].second.size() == 1)
        return radius[neighbors[completedSide].second[0]];
      else
      {
        double r = 0.0;
        using std::pow;
        for (int i = 0; i < neighbors[completedSide].second.size(); ++i)
          r += pow(radius[neighbors[completedSide].second[i]], exponent);
        r = pow(r, 1.0/exponent);
        return r;
      }
    }();

    // then check if we can proceed to the neighbor
    if (neighbors[1-completedSide].first == 1)
      for (int i = 0; i < neighbors[1-completedSide].second.size(); ++i)
        if (!visited[neighbors[1-completedSide].second[i]] && neighbors[1-completedSide].second[i] != 0)
          elemStack.push(gridView.grid().entity(seeds[neighbors[1-completedSide].second[i]]));
  };

  for (const auto& element : elements(gridView))
    for (const auto& intersection : intersections(gridView, element))
      if (intersection.boundary())
        elemStack.push(element);

  while (!elemStack.empty())
  {
    auto element = elemStack.top();
    elemStack.pop();
    handleElement(element);
  }
}

} // end namespace TreeGrowth

int main (int argc, char *argv[]) try
{
  // maybe initialize mpi
  Dune::MPIHelper::instance(argc, argv);

  // read parameter tree
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree("grid-growth.ini", params);

  TreeGrowth::Attractor<3> attractor(params.get<double>("NumAttractors", 1000), params.get<double>("MinHeight", 0.1));
  TreeGrowth::Tree<3> tree(params.get<double>("Dx", 0.01));
  Dune::VTKSequenceWriter<TreeGrowth::Tree<3>::Grid::LeafGridView> vtkWriter(tree.gridView(), "tree", "", "");
  const auto minRadius = params.get<double>("MinRadius", 0.001);
  const auto tubeLawExponent = params.get<double>("TubeLawExponent", 2.3);
  std::vector<double> radius(tree.size(), minRadius);
  vtkWriter.addCellData(radius, "radius");
  std::size_t step = 0;

  const std::size_t maxIterations = params.get<std::size_t>("MaxIterations", 1000);
  const std::size_t maxNumSegments = params.get<std::size_t>("MaxSegments", 10000);
  const std::size_t outputInterval = params.get<std::size_t>("OutputInterval", 5);
  std::size_t numSegments = tree.size();
  std::size_t iteration = 0;
  while (numSegments < maxNumSegments && iteration < maxIterations && !attractor.empty())
  {
    TreeGrowth::grow(tree, attractor, params);

    if (iteration % outputInterval == 0)
    {
      // attractor.dumpToFile("points-" + std::to_string(step) + ".txt");
      TreeGrowth::computeRadii(tree.gridView(), radius, {0, 1}, minRadius, tubeLawExponent);
      vtkWriter.write(step);
      ++step;
    }

    ++iteration;
    numSegments = tree.size();
  }

  return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
