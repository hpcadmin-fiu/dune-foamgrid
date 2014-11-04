// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONS_HH
#define DUNE_FOAMGRID_INTERSECTIONS_HH

/** \file
* \brief The FoamGridLeafIntersection and FoamGridLevelIntersection classes
*/
#include <bitset>

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

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;

    friend class FoamGridLevelIntersectionIterator<GridImp>;
    friend class FoamGridLeafIntersectionIterator<GridImp>;

    public:


        enum {dim=GridImp::dimension};

        enum {dimworld=GridImp::dimensionworld};

        typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
        typedef typename GridImp::template Codim<0>::Entity Entity;

    /**
     * \brief Initalizes an intersection.
     *
     * After initialization this object always represents the first intersection
     * related to an edge.
     * \param edge The index of the edge this intersection lives on.
     */
    FoamGridIntersection(const FoamGridEntityImp<2,dimworld>* center,
                         int edge)
        : center_(center), edgeIndex_(edge)
    {}

        //! return EntityPointer to the Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        EntityPointer inside() const {
            return FoamGridEntityPointer<0,GridImp> (center_);
        }


        //! return EntityPointer to the Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        EntityPointer outside() const {
            // Return the 'other' element on the current edge
            return FoamGridEntityPointer<0,GridImp> ((*neighbor_));
        }


        /** \brief return true if intersection is with boundary.
        */
        bool boundary () const {
           return center_->edges_[edgeIndex_]->elements_.size()==1;
        }

        //! return information about the Boundary
        int boundarySegmentIndex () const 
        {
            return center_->edges_[edgeIndex_]->boundarySegmentIndex();
        }

        //! Geometry type of an intersection
        GeometryType type () const 
        {
            return GeometryType(GeometryType::simplex, dim-1);
        }

        //! local number of codim 1 entity in self where intersection is contained in
        int indexInInside () const 
        {
            return edgeIndex_;
        }

        virtual int indexInOutside() const=0;

        //! return outer normal
        FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const 
        {
            // The intersection normal is a vector that is orthogonal to the element normal
            // and to the intersection itself.

            // only works for triangles
            assert(center_->type().isTriangle());

            // Compute vertices
            const Dune::ReferenceElement<double,dim>& refElement
                = Dune::ReferenceElements<double, dim>::general(center_->type());

            // edge vertices, oriented
            int v0 = refElement.subEntity(edgeIndex_, 1, 0, dim);
            int v1 = refElement.subEntity(edgeIndex_, 1, 1, dim);

            // opposite vertex
            int v2 = (v1+1)%3;
            if (v2==v0)
              v2 = (v0+1)%3;
            assert(v2!=v0 and v2!=v1);

            FieldVector<ctype, dimworld> edge = center_->vertex_[v0]->pos_ - center_->vertex_[v1]->pos_;
            FieldVector<ctype, dimworld> otherEdge = center_->vertex_[v2]->pos_ - center_->vertex_[v1]->pos_;

            //Cross product of edge and otherEdge is a scaled element normal
            FieldVector<ctype, dimworld> scaledElementNormal;

            if(dimworld == 3) //dimworld==3
            {
                scaledElementNormal[0] = edge[1]*otherEdge[2] - edge[2]*otherEdge[1];
                scaledElementNormal[1] = edge[2]*otherEdge[0] - edge[0]*otherEdge[2];
                scaledElementNormal[2] = edge[0]*otherEdge[1] - edge[1]*otherEdge[0];
                outerNormal_[0] = edge[1]*scaledElementNormal[2] - edge[2]*scaledElementNormal[1];
                outerNormal_[1] = edge[2]*scaledElementNormal[0] - edge[0]*scaledElementNormal[2];
                outerNormal_[2] = edge[0]*scaledElementNormal[1] - edge[1]*scaledElementNormal[0];
            } 
            else //dimworld==2
            {
                outerNormal_[0] = edge[1];
                outerNormal_[1] = -edge[0];
            }

            //Check if scaled EdgeNormal is inner normal, if yes flip
            otherEdge = center_->vertex_[v0]->pos_ - center_->vertex_[v2]->pos_;
            if(otherEdge*outerNormal_ < 0)
                outerNormal_ *= -1.0;
                        
            //Scale to unit normal vector
            //outerNormal_ /= outerNormal_.two_norm();

            return outerNormal_;
        }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const 
        {

            const Dune::ReferenceElement<double,dim>& refElement
                = Dune::ReferenceElements<double, dim>::general(center_->type());

            // edge vertices
            int v0 = refElement.subEntity(edgeIndex_, 1, 0, dim);
            int v1 = refElement.subEntity(edgeIndex_, 1, 1, dim);

            //edge length
            ctype edgeLength = (center_->vertex_[v0]->pos_- center_->vertex_[v1]->pos_).two_norm();

            FieldVector<ctype, dimworld> integrationOuterNormal_ = unitOuterNormal(local);
            integrationOuterNormal_ *= edgeLength;
           
            return integrationOuterNormal_;
        }

        //! return unit outer normal
        FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const 
        {
            unitOuterNormal_ = this->outerNormal(local);
            unitOuterNormal_ /= unitOuterNormal_.two_norm();
            return unitOuterNormal_;
        }

        //! return unit outer normal at the intersection center
        FieldVector<ctype, dimworld> centerUnitOuterNormal () const 
        {
            return unitOuterNormal(FieldVector<ctype,1>(0.5));
        }
    private:

    //! vector storing the outer normal
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;
    mutable FieldVector<typename GridImp::ctype, dimworld> unitOuterNormal_;
    mutable FieldVector<typename GridImp::ctype, dimworld> integrationOuterNormal_;

    protected:

    const FoamGridEntityImp<2,dimworld>* center_;

    /** \brief Count on which edge we are lookin' at.  */
    int edgeIndex_;

    /** \brief Iterator to the other neighbor of the intersection. */
    typename std::vector<const FoamGridEntityImp<2,dimworld>*>::const_iterator neighbor_;

};

template<class GridImp>
class FoamGridLevelIntersection
    : public FoamGridIntersection<GridImp>
{
    friend class FoamGridLevelIntersectionIterator<GridImp>;

    enum {dim=GridImp::dimension};

    // Geometry is a CachedMultiLinearGeometry
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim<1>::LocalGeometryImpl LocalGeometryImpl;

public:
    enum{ dimworld = GridImp::dimensionworld };

    FoamGridLevelIntersection(const FoamGridEntityImp<2,dimworld>* center, int edge)
                              : FoamGridIntersection<GridImp>(center, edge)
    {}

    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid level intersections are always conforming
        return true;
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
        //Not necessary 2 anymore for foamgrid t-junctions
    	//assert(this->center_->edges_[this->edgeIndex_]->elements_.size()==2);
        assert(this->neighborIndex_!=this->center_->edges_[this->edgeIndex_]->elements_.size());

        return std::find((*this->neighbor_)->edges_.begin(), (*this->neighbor_)->edges_.end(),
                         this->center_->edges_[this->edgeIndex_])
            - (*this->neighbor_)->edges_.begin();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return this->neighborIndex_!=this->center_->edges_[this->edgeIndex_]->elements_.size();
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const 
    {
        std::vector<FieldVector<double, dim> > coordinates(2);

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double, dim>::general(this->center_->type());

        coordinates[0] = refElement.position(refElement.subEntity(this->edgeIndex_, 1, 0, dim), dim);
        coordinates[1] = refElement.position(refElement.subEntity(this->edgeIndex_, 1, 1, dim), dim);
        
        geometryInInside_ = make_shared<LocalGeometryImpl>(this->type(), coordinates);

      return LocalGeometry(*geometryInInside_);
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    //! In the LevelIntersection we know that the intersection is conforming
    LocalGeometry geometryInOutside () const {
  
        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dim>& refElement
           = Dune::ReferenceElements<double, dim>::general(this->center_->type());

        std::array<FoamGridEntityImp<0,dimworld>*, 2> vtx;

        vtx[0] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 0, dim)];
        vtx[1] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 1, dim)];

        std::vector<FieldVector<double, dim> > coordinates(2);

        // Find the intersection vertices in local numbering of the outside element
        // That way we get the local orientation correctly.
        const Dune::ReferenceElement<double,dim>& refElementOther
            = Dune::ReferenceElements<double, dim>::general((*this->neighbor_)->type());

        for (int j=0; j<2; j++)
          for (int i=0; i<refElementOther.size(dim); i++)
             if (vtx[j] == (*this->neighbor_)->vertex_[refElementOther.subEntity(0,0,i,dim)])
              coordinates[j] = refElement.position(refElement.subEntity(0,0, i, dim),dim);

        geometryInOutside_ = make_shared<LocalGeometryImpl>(this->type(), coordinates);

        return LocalGeometry(*geometryInOutside_);
    }

  //! intersection of codimension 1 of this neighbor with element where iteration started.
  //! Here returned element is in GLOBAL coordinates of the element where iteration started.
  Geometry geometry () const {

        std::vector<FieldVector<double, dimworld> > coordinates(2);

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double, dim>::general(this->center_->type());

        coordinates[0] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 0, dim)]->pos_;
        coordinates[1] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 1, dim)]->pos_;
    
        geometry_ = make_shared<GeometryImpl>(this->type(), coordinates);

        return Geometry(*geometry_);
  }


    private:
    int neighborIndex_;
    //! pointer to global and local intersection geometries
    mutable Dune::shared_ptr<GeometryImpl> geometry_;
    mutable Dune::shared_ptr<LocalGeometryImpl> geometryInInside_;
    mutable Dune::shared_ptr<LocalGeometryImpl> geometryInOutside_;

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

    friend class FoamGridLeafIntersectionIterator<GridImp>;
public:

    enum {dimworld=GridImp::dimensionworld};

    enum {dim=GridImp::dimension};

    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim<1>::LocalGeometryImpl LocalGeometryImpl;

    FoamGridLeafIntersection(const FoamGridEntityImp<2,FoamGridIntersection<GridImp>::dimworld>* center,
                             int edge)
        : FoamGridIntersection<GridImp>(center, edge)
    {}

    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid leaf? intersections are always? conforming
        return true;
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
        assert(this->neighbor_!=neighborEnd_);
        // Move to the father of the edge until its level is the same as
        // the level of the neighbor
        FoamGridEntityImp<1,dimworld>* edge=(this->center_->edges_[this->edgeIndex_]);

        while(edge->level()<(*this->neighbor_)->level())
        {
            assert(edge->father_!=nullptr);
            edge=edge->father_;
        }
        assert(edge->level()==(*this->neighbor_)->level());
        assert(edge->elements_.size()==2);
        return std::find((*this->neighbor_)->edges_.begin(), (*this->neighbor_)->edges_.end(), edge)
            - (*this->neighbor_)->edges_.begin();

    }

    Geometry geometry () const 
    {
        std::vector<FieldVector<double, dimworld> > coordinates(2);

        // Get indices of two vertices of the intersection
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double, dim>::general(this->center_->type());

        // Get global coordinates of the vertices
        coordinates[0] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 0, dim)]->pos_;
        coordinates[1] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 1, dim)]->pos_;

        geometry_ = make_shared<GeometryImpl>(this->type(), coordinates);

        return Geometry(*geometry_);
      
        //DUNE_THROW(Dune::NotImplemented, "FoamGridLeafIntersection::geometry()");
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const 
    {
  
        std::vector<FieldVector<double, dim> > coordinates(2);

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double, dim>::general(this->center_->type());

        coordinates[0] = refElement.position(refElement.subEntity(this->edgeIndex_, 1, 0, dim), dim);
        coordinates[1] = refElement.position(refElement.subEntity(this->edgeIndex_, 1, 1, dim), dim);
        
        geometryInInside_ = make_shared<LocalGeometryImpl>(this->type(), coordinates);

        return LocalGeometry(*geometryInInside_);
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInOutside () const 
    {

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dim>& refElement
           = Dune::ReferenceElements<double, dim>::general(this->center_->type());

        std::array<FoamGridEntityImp<0,dimworld>*, 2> vtx;

        vtx[0] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 0, dim)];
        vtx[1] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 1, dim)];

        std::vector<FieldVector<double, dim> > coordinates(2);

        // Find the intersection vertices in local numbering of the outside element
        // That way we get the local orientation correctly.
        const Dune::ReferenceElement<double,dim>& refElementOther
            = Dune::ReferenceElements<double, dim>::general((*this->neighbor_)->type());
          
        for (int j=0; j<2; j++) //check both vertices of the considered edge
         for (int i=0; i<refElementOther.size(dim); i++) //against all vertices of neighbor element
            if (vtx[j] == (*this->neighbor_)->vertex_[refElementOther.subEntity(0, 0, i, dim)])
              coordinates[j] = refElement.position(refElement.subEntity(0, 0, i, dim), dim);

        geometryInOutside_ = make_shared<LocalGeometryImpl>(this->type(), coordinates);

        return LocalGeometry(*geometryInOutside_);
    }

    //! return true if across the edge a neighbor on this level exists
    bool neighbor () const {
        return this->neighbor_!=neighborEnd_;
    }

private:

    FoamGridEntityImp<1,dimworld>* edgePointer_;
    /** \brief Iterator to the other neighbor of the intersection. */
    typename std::vector<const FoamGridEntityImp<2,dimworld>*>::const_iterator neighborEnd_;
    //! pointer to global and local intersection geometries
    mutable Dune::shared_ptr<GeometryImpl> geometry_;
    mutable Dune::shared_ptr<LocalGeometryImpl> geometryInInside_;
    mutable Dune::shared_ptr<LocalGeometryImpl> geometryInOutside_;
};


}  // namespace Dune

#endif
