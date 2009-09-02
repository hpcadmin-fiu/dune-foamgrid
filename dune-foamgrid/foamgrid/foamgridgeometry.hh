#ifndef DUNE_IDENTITYGRIDGEOMETRY_HH
#define DUNE_IDENTITYGRIDGEOMETRY_HH

/** \file
* \brief The FoamGridGeometry class and its specializations
*/

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>

namespace Dune {

template<int mydim, int coorddim, class GridImp>
class FoamGridGeometry :
        public GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> >
{

    typedef typename GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> > Base;

    public:

    /** \brief Constructor with a geometry type and a set of corners */
    void setup(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates)
    {
        // set up base class
        // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
        Base::operator=(Base(type,coordinates));
    }

};


}  // namespace Dune

#endif
