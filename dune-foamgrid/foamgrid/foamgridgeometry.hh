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
        //public GeometryDefaultImplementation <mydim, coorddim, GridImp, FoamGridGeometry>
        public GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> >
{

    typedef typename GenericGeometry::BasicGeometry<mydim, GenericGeometry::DefaultGeometryTraits<typename GridImp::ctype,mydim,coorddim> > Base;

    public:

    /** \brief Constructor with a geometry type and a set of corners */
    void setup(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,2> >& coordinates)
    {
        // set up base class
        // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
        Base::operator=(Base(type,coordinates));
    }

#if 0
    private:
    
        typedef typename GridImp::ctype ctype;
    
    
    public:
    
        // The codimension of this entitypointer wrt the host grid
        enum {CodimInHostGrid = GridImp::HostGridType::dimension - mydim};
        enum {DimensionWorld = GridImp::HostGridType::dimensionworld};
    
        // select appropiate hostgrid geometry via typeswitch
        typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridGeometryType;
        typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridLocalGeometryType;
        
        typedef typename SelectType<coorddim==DimensionWorld, HostGridGeometryType, HostGridLocalGeometryType>::Type HostGridGeometry;
        
        
        /** Default constructor.
        */
        FoamGridGeometry(const HostGridGeometry& hostGeometry)
            : hostGeometry_(hostGeometry)
        {
        }
    
        
        /** \brief Return the element type identifier
        */
        GeometryType type () const {
            return hostGeometry_.type();
        }
    
        
        //! return the number of corners of this element. Corners are numbered 0...n-1
        int corners () const {
            return hostGeometry_.corners();
        }
    
        
        //! access to coordinates of corners. Index is the number of the corner
        const FieldVector<ctype, coorddim> corner (int i) const {
            return hostGeometry_.corner(i);
        }
    
        
        /** \brief Maps a local coordinate within reference element to
        * global coordinate in element  */
        FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const{
            return hostGeometry_.global(local);
        }
    
        
        /** \brief Maps a global coordinate within the element to a
        * local coordinate in its reference element */
        FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const {
            return hostGeometry_.local(global);
        }
    
        
        //! Returns true if the point is in the current element
        bool checkInside(const FieldVector<ctype, mydim> &local) const {
            return hostGeometry_.checkInside(local);
        }
    
        
        /**
        */
        ctype integrationElement (const FieldVector<ctype, mydim>& local) const {
            return hostGeometry_.integrationElement(local);
        }
    
        
        //! The Jacobian matrix of the mapping from the reference element to this element
        const FieldMatrix<ctype, coorddim,mydim>& jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
            return hostGeometry_.jacobianInverseTransposed(local);
        }
    
        
        const HostGridGeometry& hostGeometry_;
#endif
};


}  // namespace Dune

#endif
