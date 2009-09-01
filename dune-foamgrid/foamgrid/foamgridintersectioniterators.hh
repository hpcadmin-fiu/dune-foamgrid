#ifndef DUNE_FOAMGRID_INTERSECTIONITERATORS_HH
#define DUNE_FOAMGRID_INTERSECTIONITERATORS_HH

#include <dune-foamgrid/foamgrid/foamgridintersections.hh>

/** \file
* \brief The FoamGridLeafIntersectionIterator and FoamGridLevelIntersectionIterator classes
*/

namespace Dune {

/** \brief Iterator over all element neighbors
* \ingroup FoamGrid
* Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
* a neighbor is an entity of codimension 0 which has a common entity of codimension 1
* These neighbors are accessed via a IntersectionIterator. This allows the implement
* non-matching meshes. The number of neighbors may be different from the number
* of an element!
*/
template<class GridImp>
class FoamGridLeafIntersectionIterator
{
    
    enum {dim=GridImp::dimension};
    
    enum {dimworld=GridImp::dimensionworld};
    
    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;
    
    typedef typename GridImp::HostGridType::template Codim<0>::Entity::LeafIntersectionIterator HostLeafIntersectionIterator;
    
public:
    
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef Dune::Intersection<const GridImp, Dune::FoamGridLeafIntersectionIterator> Intersection;
    
    FoamGridLeafIntersectionIterator(const GridImp* identityGrid,
                                         const HostLeafIntersectionIterator& hostIterator)
        : selfLocal_(NULL), neighborLocal_(NULL), intersectionGlobal_(NULL),
          identityGrid_(identityGrid), 
          hostIterator_(hostIterator)
    {}
        
    //! The Destructor
    ~FoamGridLeafIntersectionIterator() {};
    
    //! equality
    bool equals(const FoamGridLeafIntersectionIterator<GridImp>& other) const {
        return hostIterator_ == other.hostIterator_;
    }

    
    //! prefix increment
    void increment() {
        ++hostIterator_;

        // Delete intersection geometry objects, if present
        if (intersectionGlobal_ != NULL) {
            delete intersectionGlobal_;
            intersectionGlobal_ = NULL;
        }
        
        if (selfLocal_ != NULL) {
            delete selfLocal_;
            selfLocal_ = NULL;
        }
        
        if (neighborLocal_ != NULL) {
            delete neighborLocal_;
            neighborLocal_ = NULL;
        }
    }

    //! \brief dereferencing
    const Intersection & dereference() const {
        return reinterpret_cast<const Intersection&>(*this);
    }
    
private:
        //**********************************************************
        //  private methods
        //**********************************************************

    //! pointer to element holding the selfLocal and selfGlobal information.
    //! This element is created on demand.
    mutable MakeableInterfaceObject<LocalGeometry>* selfLocal_;
    mutable MakeableInterfaceObject<LocalGeometry>* neighborLocal_;
    
    //! pointer to element holding the neighbor_global and neighbor_local
    //! information.
    mutable MakeableInterfaceObject<Geometry>* intersectionGlobal_;

    const GridImp* identityGrid_;

    HostLeafIntersectionIterator hostIterator_;
};




//! \todo Please doc me !
template<class GridImp>
class FoamGridLevelIntersectionIterator
{
    
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };

    //friend class OneDGridEntity<0,dim,GridImp>;

    //! Constructor for a given grid entity and a given neighbor
    FoamGridLevelIntersectionIterator(FoamGridElement* center, int nb) 
        : intersection_(center,nb)
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLevelIntersectionIterator(FoamGridElement* center) 
        : intersection_(center,center->corners())
    {}

public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef Dune::Intersection<const GridImp, Dune::FoamGridLevelIntersection> Intersection;

  //! equality
  bool equals(const FoamGridLevelIntersectionIterator<GridImp>& other) const {
      return (center_ == other.center_) && (neighbor_ == other.neighbor_);
  }

    //! prefix increment
    void increment() {
        neighbor_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
        return reinterpret_cast<const Intersection&>(*this);
    }

    FoamGridElement* target() const {
        const bool isValid = center_ && neighbor_>=0 && neighbor_<2;

        if (!isValid)
            return center_;
        else if (neighbor_==0) 
            return center_->pred_;
        else 
            return center_->succ_;

    }

private:
  //**********************************************************
  //  private methods 
  //**********************************************************

    /** \brief The actual intersection
    */
    mutable MakeableInterfaceObject<FoamGridLevelIntersection<GridImp> > intersection_;

};


}  // namespace Dune

#endif
