// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONS_HH
#define DUNE_FOAMGRID_INTERSECTIONS_HH

/** \file
* \brief The FoamGridLeafIntersection and FoamGridLevelIntersection classes
*/

#include <dune/grid/common/intersection.hh>

#include <dune/foamgrid/foamgrid/foamgridintersectioniterators.hh>
#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridgeometry.hh>
#include <dune/foamgrid/foamgrid/foamgridentitypointer.hh>

namespace Dune {

template <class GridImp>
class FoamGridLevelIntersectionIterator;


//! \brief Base class of all intersections within FoamGrid
//!
//! encapsulates common functionality of level and leaf intersections.
template<class GridImp>
class FoamGridIntersection
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

    /**
     * \brief Initalizes an intersection.
     *
     * After initialization this object always represents the first intersection 
     * related to an edge.
     * \param edge The index of the edge this intersection lives on.
     */
    FoamGridIntersection(const FoamGridEntityImp<2,dimworld>* center, int edge)
        : center_(center), edgeIndex_(edge), neighbor_(0)
    {
        if(edge!=center_->corners() &&
           center_->edges_[edgeIndex_]->elements_.size() &&
           center_==center_->edges_[edgeIndex_]->elements_[0])
            // Move index to point to the first real neighbor 
            ++neighbor_;
    }

    void increment(){
        ++neighbor_;
        while(neighbor_<center_->edges_[edgeIndex_]->elements_.size() &&
              center_==center_->edges_[edgeIndex_]->elements_[neighbor_])
            ++neighbor_;
        
        if(neighbor_<center_->edges_[edgeIndex_]->elements_.size())
            return;
        ++edgeIndex_;
        neighbor_=0;
        if(edgeIndex_!=center_->corners() &&
           center_==center_->edges_[edgeIndex_]->elements_[0])
            // Move index to point to the first real neighbor 
            ++neighbor_;
    }
    
        //! return EntityPointer to the Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        EntityPointer inside() const {
            return FoamGridEntityPointer<0,GridImp> (center_);
        }

        
        //! return EntityPointer to the Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        EntityPointer outside() const {
            // Return the 'other' element on the current edge
            return FoamGridEntityPointer<0,GridImp> (center_->edges_[edgeIndex_]->elements_[neighbor_]);
        }
        
        
    /** \brief return true if intersection is with boundary.
    */
    bool boundary () const {
        return center_->edges_[edgeIndex_]->elements_.size()==1;
    }
        
        
    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
        return center_->edges_[edgeIndex_]->elements_.size()>1;
    }
    
    
    //! return information about the Boundary
    int boundaryId () const DUNE_DEPRECATED {
        return center_->edges_[edgeIndex_]->boundarySegmentIndex();
    }
        
    //! return information about the Boundary
    int boundarySegmentIndex () const {
        return center_->edges_[edgeIndex_]->boundarySegmentIndex();
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
        LocalGeometry geometryInInside () const {

            std::vector<FieldVector<double, dim> > coordinates(2);

            // Get two vertices of the intersection
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            coordinates[0] = refElement.position(refElement.subEntity(edgeIndex_, 1, 0, dim),dim);
            coordinates[1] = refElement.position(refElement.subEntity(edgeIndex_, 1, 1, dim),dim);

            return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(type(), coordinates));
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in LOCAL coordinates of neighbor
        LocalGeometry geometryInOutside () const {

            std::vector<FieldVector<double, dim> > coordinates(2);

            // Get two vertices of the intersection
            const FoamGridEntityImp<2,dimworld>* outside = center_->edges_[edgeIndex_]->elements_[neighbor_];
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(outside->type());

            int idxInOutside = indexInOutside();
            
            coordinates[0] = refElement.position(refElement.subEntity(idxInOutside, 1, 0, dim),dim);
            coordinates[1] = refElement.position(refElement.subEntity(idxInOutside, 1, 1, dim),dim);

            return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(type(), coordinates));
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in GLOBAL coordinates of the element where iteration started.
        Geometry geometry () const {

            std::vector<FieldVector<double, dimworld> > coordinates(2);

            // Get two vertices of the intersection
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            coordinates[0] = center_->vertex_[refElement.subEntity(edgeIndex_, 1, 0, dim)]->pos_;
            coordinates[1] = center_->vertex_[refElement.subEntity(edgeIndex_, 1, 1, dim)]->pos_;

            return Geometry(FoamGridGeometry<dim-1, dimworld, GridImp>(type(), coordinates));
        }
        
        
        //! local number of codim 1 entity in self where intersection is contained in
        int indexInInside () const {
            return edgeIndex_;
        }
        
        
        //! local number of codim 1 entity in neighbor where intersection is contained
        int indexInOutside () const {
            assert(center_->edges_[edgeIndex_]->elements_.size()==2);
            const FoamGridEntityImp<2,dimworld>* other = center_->edges_[edgeIndex_]->elements_[neighbor_];
            assert(other);
            
            return std::find(other->edges_.begin(), other->edges_.end(), center_->edges_[edgeIndex_]) 
                - other->edges_.begin();
        }
        
          
        //! return outer normal
        FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {
            // The intersection normal is a vector that is orthogonal to the element normal
            // and to the intersection itself.
            
            // only works for triangles
            assert(center_->type().isTriangle());

            // Compute vertices
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            // edge vertices, oriented
            int v0 = std::min(refElement.subEntity(edgeIndex_, 1, 0, dim), refElement.subEntity(edgeIndex_, 1, 1, dim));
            int v1 = std::max(refElement.subEntity(edgeIndex_, 1, 0, dim), refElement.subEntity(edgeIndex_, 1, 1, dim));

            // opposite vertex
            int v2 = (v1+1)%3;

            // Compute oriented edge
            FieldVector<ctype, dimworld> edge = center_->vertex_[v1]->pos_ - center_->vertex_[v0]->pos_;

            // compute triangle edge normal
            FieldVector<ctype, dimworld> scaledEdge = edge;
            edge *= edge*(center_->vertex_[v2]->pos_ - center_->vertex_[v0]->pos_);
            FieldVector<ctype, dimworld> normal = center_->vertex_[v2]->pos_ - center_->vertex_[v0]->pos_;
            normal -= scaledEdge;
            normal *= -1;
            return normal;
        }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {

            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            ctype edgeLength = (center_->vertex_[refElement.subEntity(edgeIndex_, 1, 0, dim)]->pos_
                                - center_->vertex_[refElement.subEntity(edgeIndex_, 1, 1, dim)]->pos_).two_norm();

            FieldVector<ctype, dimworld> normal = unitOuterNormal(local);
            normal *= edgeLength;
            return normal;
        }

        //! return unit outer normal
        FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {
            FieldVector<ctype, dimworld> outerNormal = this->outerNormal(local);
            outerNormal /= outerNormal.two_norm();
            return outerNormal;
        }

        //! return unit outer normal at the intersection center
        FieldVector<ctype, dimworld> centerUnitOuterNormal () const {
            FieldVector<ctype, dimworld> outerNormal = this->outerNormal(FieldVector<ctype,1>(0.5));
            outerNormal /= outerNormal.two_norm();
            return outerNormal;
        }

    private:

    const FoamGridEntityImp<2,dimworld>* center_;
 
    //! vector storing the outer normal 
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;

    /** \brief Count on which edge we are lookin' at.  */
    int edgeIndex_;
    
    /** 
     * \brief Count on which neighbor of the edge we are looking at.
     *
     * At a t-junction different intersections might occur for the same edge.
     * In this case each intersection has a different set of neighbor elements.
     * This index is always the one of the current neighbor in FoamGridEdge::elements_
     */
    int neighbor_;
    
};

template<class GridImp>
class FoamGridLevelIntersection
    : public FoamGridIntersection<GridImp>
{

public:
    enum{ dimworld = GridImp::dimensionworld };
    
    FoamGridLevelIntersection(const FoamGridEntityImp<2,dimworld>* center, int edge)
        : FoamGridIntersection<GridImp>(center, edge)
    {}
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
    : public FoamGridIntersection<GridImp>
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
        : selfLocal_(nullptr), neighborLocal_(nullptr), intersectionGlobal_(nullptr),
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
        int boundarySegmentIndex () const {
            return hostIterator_->boundarySegmentIndex();
        }

        //! return information about the Boundary
        int boundaryId () const DUNE_DEPRECATED {
            return hostIterator_->boundarySegmentIndex();
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
            if (selfLocal_ == nullptr)
                selfLocal_ = new MakeableInterfaceObject<LocalGeometry>(hostIterator_->intersectionSelfLocal());
                
            return *selfLocal_;
        }
    
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in LOCAL coordinates of neighbor
        const  LocalGeometry& geometryInOutside () const {
            if (neighborLocal_ == nullptr)
                neighborLocal_ = new MakeableInterfaceObject<LocalGeometry>(hostIterator_->intersectionNeighborLocal());
                
            return *neighborLocal_;
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in GLOBAL coordinates of the element where iteration started.
        const  Geometry& geometry () const {
            if (intersectionGlobal_ == nullptr)
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
        
    //! return unit outer normal at the intersection center
    FieldVector<ctype, dimworld> centerUnitOuterNormal () const {
        FieldVector<ctype, dimworld> outerNormal = this->outerNormal(FieldVector<ctype,1>(0.5));
        outerNormal /= outerNormal.two_norm();
        return outerNormal;
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
