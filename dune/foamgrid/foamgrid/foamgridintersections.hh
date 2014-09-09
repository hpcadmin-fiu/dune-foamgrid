// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONS_HH
#define DUNE_FOAMGRID_INTERSECTIONS_HH

/** \file
* \brief The FoamGridLeafIntersection and FoamGridLevelIntersection classes
*/

#include <dune/grid/common/intersection.hh>

#include "foamgridintersectioniterators.hh"
#include "foamgridvertex.hh"
#include "foamgridgeometry.hh"
#include "foamgridentitypointer.hh"

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

    FoamGridIntersection(const FoamGridEntityImp<1,dimworld>* center,
				  int vertex)
        : center_(center), vertexIndex_(vertex)
    {
    }

        //! return EntityPointer to the Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        EntityPointer inside() const {
            return FoamGridEntityPointer<0,GridImp> (center_);
        }


        //! return EntityPointer to the Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        EntityPointer outside(int id=0) const {
            // Return the 'other' element on the current vertex
            return FoamGridEntityPointer<0,GridImp> (*(neighbors_[id]));
        }


    /** \brief return true if intersection is with boundary.
    */
    bool boundary () const {
	return center_->vertex_[vertexIndex_]->elements_.size()==1;
    }

    /** \brief return the number of neighbors
     */
    int NumOfNeighbors () const {
      return center_->vertex_[vertexIndex_]->elements_.size();
    }

    //! return information about the Boundary
    int boundarySegmentIndex () const {
        return center_->vertex_[vertexIndex_]->boundarySegmentIndex();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
        return GeometryType(GeometryType::simplex, dim-1);
    }

         //! local number of the intersection in the edge that was given in the constructor
        int indexInInside () const {
            return vertexIndex_;
        }

        virtual int indexInOutside() const=0;

        //! return outer normal
        FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {

	    FieldVector<ctype, dimworld> edge = center_->vertex_[vertexIndex_]->pos - center_->vertexIndex_[1-vertexIndex_]->pos;
	    FieldVector<ctype, dimworld> normal = edge /= edge.two_norm;

            return normal;
       }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {

            return this->outerNormal(local);
        }

        //! return unit outer normal
        FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {

            return this->outerNormal(local);
        }

        //! return unit outer normal at the intersection center
        FieldVector<ctype, dimworld> centerUnitOuterNormal () const {

            return this->unitOuterNormal(FieldVector<ctype,0>(0.5));

        }
    private:

    //! vector storing the outer normal
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;

    protected:

    const FoamGridEntityImp<1,dimworld>* center_;

    /** \brief Count on which vertex we are lookin' at.  */

    int vertexIndex_;

    /** \brief Iterator to the neighbors of the intersection. */
    std::vector<typename std::vector<const FoamGridEntityImp<1,dimworld>*>::const_iterator > neighbors_;

};

template<class GridImp>
class FoamGridLevelIntersection
    : public FoamGridIntersection<GridImp>
{
    friend class FoamGridLevelIntersectionIterator<GridImp>;

    enum {dim=GridImp::dimension};

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

public:
    enum{ dimworld = GridImp::dimensionworld };

    FoamGridLevelIntersection(const FoamGridEntityImp<1,dimworld>* center, int vertex)
                              : FoamGridIntersection<GridImp>(center, vertex)
    {}

    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid level intersections are always conforming
        return true;
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
//     int indexInOutside (int id=0) const {
//         assert(this->center_->vertex_[this->vertexIndex_]->elements_.size()>=2);
//         assert(this->neighborIndex_!=this->center_->vertex_[this->vertexIndex_]->elements_.size());
//
//         return std::find((*this->neighbors_[id])->vertex_.begin(), (*this->neighbors_[id])->vertex_.end(),
//                          this->center_->vertex_[this->vertexIndex_])
//             - (*this->neighbors_[id])->vertex_.begin();
//     }

    //! return true if across the edge an neighbor on this level exists
//     bool neighbor () const {
//       return this->neighborIndex_!=this->center_->vertex_[this->vertexIndex_]->elements_.size();
//     }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
//     LocalGeometry geometryInInside () const {
//
//       std::vector<FieldVector<double, dim> > coordinates(1);
//
//       coordinates[0] = this->VertexIndex_;
//
//
//       return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(this->type(), coordinates));
//     }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    //! In the LevelIntersection we know that the intersection is conforming
//     LocalGeometry geometryInOutside () const {
//
//     // Get two vertices of the intersection
//     const Dune::ReferenceElement<double,dim>& refElement
//         = Dune::ReferenceElements<double, dim>::general(this->center_->type());
//
//     std::array<FoamGridEntityImp<0,dimworld>*, 2> vtx;
//
//     vtx[0] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 0, dim)];
//     vtx[1] = this->center_->vertex_[refElement.subEntity(this->edgeIndex_, 1, 1, dim)];
//
//     std::vector<FieldVector<double, dim> > coordinates(2);
//
//     // Find the intersection vertices in local numbering of the outside element
//     // That way we get the local orientation correctly.
//     const Dune::ReferenceElement<double,dim>& refElementOther
//         = Dune::ReferenceElements<double, dim>::general((*this->neighbor_)->type());
//
//     for (int i=0; i<refElementOther.size(dim); i++)
//       for (int j=0; j<2; j++)
//         if (vtx[j] == (*this->neighbor_)->vertex_[refElementOther.subEntity(0,0,i,dim)])
//           coordinates[j] = refElement.position(refElement.subEntity(0,0, i, dim),dim);
//
//     return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(this->type(), coordinates));
//   }

  //! intersection of codimension 1 of this neighbor with element where iteration started.
  //! Here returned element is in GLOBAL coordinates of the element where iteration started.
//   Geometry geometry () const {
//
//     std::vector<FieldVector<double, dimworld> > coordinates(1);;
//     coordinates[0] = this->center_->vertex_[this->VertexIndex_]->pos_;
//
//     return Geometry(FoamGridGeometry<dim-1, dimworld, GridImp>(this->type(), coordinates));
//   }


    private:
    int neighborIndex_;

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

    FoamGridLeafIntersection(const FoamGridEntityImp<1,FoamGridIntersection<GridImp>::dimworld>* center,
                             int vertex)
        : FoamGridIntersection<GridImp>(center, vertex)
    {}

    //! Return true if this is a conforming intersection
    bool conforming () const {
        DUNE_THROW(Dune::NotImplemented, "FoamGridLeafIntersection::conforming");
    }

    //! local number of the intersection in the neighbor
//     int indexInOutside (int id=0) const {
//         assert(this->neighbors_[id]!=neighborEnd_);
//         // Move to the father of the edge until its level is the same as
//         // the level of the neighbor
//         FoamGridEntityImp<1,dimworld>* edge=this->center_;
//
//         while(edge->level()<(*this->neighbors_[id])->level())
//         {
//             assert(edge->father_!=nullptr);
//             edge=edge->father_;
//         }
//
//         assert(edge->level()==(*this->neighbors_[id])->level());
//         assert(edge->elements_.size()>=2);
//         return std::find((*this->neighbors_[id])->vertex_.begin(), (*this->neighbors_[id])->vertex_.end(), edge)
//             - (*this->neighbors_[id])->vertex_.begin();
//
//     }

    Geometry geometry () const {
      DUNE_THROW(Dune::NotImplemented, "FoamGridLeafIntersection::geometry()");
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
//     LocalGeometry geometryInInside () const {
//
//         std::vector<FieldVector<double, dim> > coordinates(1);
//
//         coordinates[0]=this->center_->vertex_[this->vertexIndex_]->pos_;
//
//         return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(this->type(), coordinates));
//     }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
//     LocalGeometry geometryInOutside () const {
//         /*
//         if(edgePointer_.level()==(*this->neighbor_)->level())
//             return FoamGridIntersection<GridImp>::geometryInOutside();
//         */
//         std::vector<FieldVector<double, dim> > coordinates(2);
//
//         coordinates[0]=(*this->neighbor_)->globalToLocal(edgePointer_->vertex_[0]->pos_);
//
//         return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(this->type(), coordinates));
//     }

     //! return outer normal multiplied by the integration element
//         FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {
//             ctype edgeLength = (edgePointer_->vertex_[0]->pos_ - edgePointer_->vertex_[1]->pos_).two_norm();
//             FieldVector<ctype, dimworld> normal = this->unitOuterNormal(local);
//             normal *= edgeLength;
//             return normal;
//         }
    //! return true if across the edge an neighbor on this level exists
    bool neighbor (int id=0) const {
        return this->neighbors_[id]!=neighborEnd_;
    }
private:
    FoamGridEntityImp<1,dimworld>* edgePointer_;
    /** \brief Iterator to the other neighbor of the intersection. */
    typename std::vector<const FoamGridEntityImp<1,dimworld>*>::const_iterator neighborEnd_;
};


}  // namespace Dune

#endif
