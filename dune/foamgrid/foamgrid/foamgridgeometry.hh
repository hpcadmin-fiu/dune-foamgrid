#ifndef DUNE_FOAMGRID_GEOMETRY_HH
#define DUNE_FOAMGRID_GEOMETRY_HH

/** \file
* \brief The FoamGridGeometry class
*/

#include <vector>

#include <dune/common/version.hh>
#include <dune/geometry/affinegeometry.hh>


namespace Dune {

#if ! DUNE_VERSION_NEWER(DUNE_GRID,2,4)
template<int mydim, int coorddim, class GridImp> class FoamGridGeometry;

namespace FacadeOptions
{
    template< int mydim, int cdim, class GridImp>
    struct StoreGeometryReference<mydim, cdim, GridImp, FoamGridGeometry>
    {
        static const bool v = false;
    };
}
#endif

template<int mydim, int coorddim, class GridImp>
class FoamGridGeometry :
        public AffineGeometry<typename GridImp::ctype, mydim, coorddim>
{
    public:

    /**
     * \brief This is DefaultConstructor
     */
    FoamGridGeometry() {}

    /**
     * \brief Construct geometry from coordinate vector
     */
    FoamGridGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates) :
        AffineGeometry<typename GridImp::ctype, mydim, coorddim>(type, coordinates)
    {}

};


}  // namespace Dune

#endif
