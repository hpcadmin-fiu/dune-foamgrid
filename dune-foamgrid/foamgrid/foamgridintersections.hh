#ifndef DUNE_FOAMGRID_INTERSECTIONS_HH
#define DUNE_FOAMGRID_INTERSECTIONS_HH

/** \file
* \brief The FoamGridLeafIntersection and FoamGridLevelIntersection classes
*/

namespace Dune {




//! \todo Please doc me !
template<class GridImp>
class FoamGridLevelIntersection
{
    
        enum {dim=GridImp::dimension};
    
        enum {dimworld=GridImp::dimensionworld};
    
        // The type used to store coordinates
        typedef typename GridImp::ctype ctype;

    friend class FoamGridLevelIntersectionIterator<GridImp>;
    
    public:

        typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
        typedef typename GridImp::template Codim<1>::Geometry Geometry;
        typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
        typedef typename GridImp::template Codim<0>::Entity Entity;
        typedef Dune::Intersection<const GridImp, Dune::FoamGridLevelIntersectionIterator> Intersection;

    FoamGridLevelIntersection(const FoamGridElement* center, int nb)
        : center_(center), neighbor_(nb),
          geometryInInside_(FoamGridGeometry<dim-1,dim,GridImp>()),
          geometryInOutside_(FoamGridGeometry<dim-1,dim,GridImp>()),
          geometry_(FoamGridGeometry<dim-1,dimworld,GridImp>())
    {}

        //! return EntityPointer to the Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        EntityPointer inside() const {
            return FoamGridEntityPointer<0,GridImp> (center_);
        }

        
        //! return EntityPointer to the Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        EntityPointer outside() const {
            assert(center_->edges_[neighbor_]->elements_.size()==2);
            // Return the 'other' element on the current edge
            return (center_->edges_[neighbor_]->elements_[0]==center_)
                ? FoamGridEntityPointer<0,GridImp> (center_->edges_[neighbor_]->elements_[1])
                : FoamGridEntityPointer<0,GridImp> (center_->edges_[neighbor_]->elements_[0]);
        }
        
        
    /** \brief return true if intersection is with boundary.
    */
    bool boundary () const {
        return center_->edges_[neighbor_]->elements_.size()==1;
    }
        
        
    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
        return center_->edges_[neighbor_]->elements_.size()>1;
    }
    
    
    //! return information about the Boundary
    int boundaryId () const {
        return center_->edges_[neighbor_]->boundaryId();
    }
        
    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid level intersections are always conforming
        return true;
    }
        
    //! Geometry type of an intersection
    GeometryType type () const {
        return GeometryType(GeometryType::simplex, dim-1);
    }


        //! intersection of codimension 1 of this neighbor with element where
        //! iteration started.
        //! Here returned element is in LOCAL coordinates of the element
        //! where iteration started.
        const LocalGeometry& geometryInInside () const {

            std::vector<FieldVector<double, dim> > coordinates(2);

            GridImp::getRealImplementation(geometryInInside_).setup(GeometryType(GeometryType::simplex,dim-1), coordinates);

            return geometryInInside_;
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in LOCAL coordinates of neighbor
        const  LocalGeometry& geometryInOutside () const {

            std::vector<FieldVector<double, dim> > coordinates(2);

            GridImp::getRealImplementation(geometryInOutside_).setup(GeometryType(GeometryType::simplex,dim-1), coordinates);
                
            return geometryInOutside_;
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in GLOBAL coordinates of the element where iteration started.
        const Geometry& geometry () const {

            std::vector<FieldVector<double, dimworld> > coordinates(2);

            GridImp::getRealImplementation(geometry_).setup(GeometryType(GeometryType::simplex,dim-1), coordinates);
                
            return geometry_;
        }
        
        
        //! local number of codim 1 entity in self where intersection is contained in
        int indexInInside () const {
            return neighbor_;
        }
        
        
        //! local number of codim 1 entity in neighbor where intersection is contained
        int indexInOutside () const {
            DUNE_THROW(NotImplemented, "!");
            return 0;
        }
        
          
        //! return outer normal
        FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {
            DUNE_THROW(NotImplemented, "!");
            return FieldVector<ctype, dimworld>(0);
        }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {
            DUNE_THROW(NotImplemented, "!");
            return FieldVector<ctype, dimworld>(0);
        }

        //! return unit outer normal
        FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {
            FieldVector<ctype, dimworld> outerNormal = outerNormal(local);
            outerNormal /= outerNormal.two_norm();
            return outerNormal;
        }

    private:

    const FoamGridElement* center_;
 
    //! vector storing the outer normal 
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;

    /** \brief Count on which neighbor we are lookin' at.  */
    int neighbor_;

    /** \brief The geometry that's being returned when intersectionSelfLocal() is called
    */
    mutable MakeableInterfaceObject<LocalGeometry> geometryInInside_;

    /** \brief The geometry that's being returned when intersectionNeighborLocal() is called
    */
    mutable MakeableInterfaceObject<LocalGeometry> geometryInOutside_;
    
    //! The geometry that's being returned when intersectionSelfGlobal() is called
    mutable MakeableInterfaceObject<Geometry> geometry_;
};



/** \brief Iterator over all element neighbors
* \ingroup FoamGrid
* Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
* a neighbor is an entity of codimension 0 which has a common entity of codimension 1
* These neighbors are accessed via a IntersectionIterator. This allows the implement
* non-matching meshes. The number of neighbors may be different from the number
* of an element!
*/
template<class GridImp>
class FoamGridLeafIntersection
    : public FoamGridLevelIntersection<GridImp>
{
#if 0
    enum {dim=GridImp::dimension};
    
    enum {dimworld=GridImp::dimensionworld};
    
    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;
    
#endif
public:
    
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef Dune::Intersection<const GridImp, Dune::FoamGridLeafIntersection> Intersection;

#if 0    
    FoamGridLeafIntersection(const GridImp* identityGrid,
                                         const HostLeafIntersectionIterator& hostIterator)
        : selfLocal_(NULL), neighborLocal_(NULL), intersectionGlobal_(NULL),
          identityGrid_(identityGrid), 
          hostIterator_(hostIterator)
    {}
        
    //! The Destructor
    ~FoamGridLeafIntersection() {};
    
    //! return EntityPointer to the Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        EntityPointer inside() const {
            return FoamGridEntityPointer<0,GridImp> (identityGrid_, hostIterator_->inside());
        }

    
        //! return EntityPointer to the Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        EntityPointer outside() const {
            return FoamGridEntityPointer<0,GridImp> (identityGrid_, hostIterator_->outside());
        }

    
        //! return true if intersection is with boundary.
        bool boundary () const {
            return hostIterator_->boundary();
        }
    
        
        //! return true if across the edge an neighbor on this level exists
        bool neighbor () const {
            return hostIterator_->neighbor();
        }
        
        
        //! return information about the Boundary
        int boundaryId () const {
            return hostIterator_->boundaryId();
        }

    //! Return true if this is a conforming intersection
    bool conforming () const {
        return hostIterator_->conforming();
    }
        
    //! Geometry type of an intersection
    GeometryType type () const {
        return hostIterator_->type();
    }


        //! intersection of codimension 1 of this neighbor with element where
        //! iteration started.
        //! Here returned element is in LOCAL coordinates of the element
        //! where iteration started.
        const  LocalGeometry& geometryInInside () const {
            if (selfLocal_ == NULL)
                selfLocal_ = new MakeableInterfaceObject<LocalGeometry>(hostIterator_->intersectionSelfLocal());
                
            return *selfLocal_;
        }
    
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in LOCAL coordinates of neighbor
        const  LocalGeometry& geometryInOutside () const {
            if (neighborLocal_ == NULL)
                neighborLocal_ = new MakeableInterfaceObject<LocalGeometry>(hostIterator_->intersectionNeighborLocal());
                
            return *neighborLocal_;
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in GLOBAL coordinates of the element where iteration started.
        const  Geometry& geometry () const {
            if (intersectionGlobal_ == NULL)
                intersectionGlobal_ = new MakeableInterfaceObject<Geometry>(hostIterator_->intersectionGlobal());
                
            return *intersectionGlobal_;
        }
    
        
        //! local number of codim 1 entity in self where intersection is contained in
        int indexInInside () const {
            return hostIterator_->indexInInside();
        }
    
        
        //! local number of codim 1 entity in neighbor where intersection is contained
        int indexInOutside () const {
            return hostIterator_->indexInOutside();
        }
    
    
        //! return outer normal
        FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
            return hostIterator_->outerNormal(local);
        }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
            return hostIterator_->integrationOuterNormal(local);
        }

        //! return unit outer normal
        FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
            return hostIterator_->unitOuterNormal(local);
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
#endif
};


}  // namespace Dune

#endif
