// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONITERATORS_HH
#define DUNE_FOAMGRID_INTERSECTIONITERATORS_HH

#include <dune/foamgrid/foamgrid/foamgridintersections.hh>
#include <dune/foamgrid/foamgrid/foamgridentity.hh>
#include <dune/foamgrid/foamgrid/foamgridvertex.hh>

#include <map>
#include <dune/common/shared_ptr.hh>

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
* of corners of an element!
*/
template<class GridImp>
class FoamGridLeafIntersectionIterator
{

    enum {dimworld = GridImp::dimensionworld};
    enum {dimgrid  = GridImp::dimension};

    typedef std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*> ElementVector;

    typedef std::map<int, ElementVector> MapType;

    // Only the codim-0 entity is allowed to call the constructors
    friend class FoamGridEntity<0,dimgrid,GridImp>;

    template<typename, typename, typename>
    friend class Dune::IntersectionIterator;

    FoamGridLeafIntersectionIterator()
    {}

    //! Constructor for a given grid entity and a given neighbor
    FoamGridLeafIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center, int facet)
        : intersection_(FoamGridLeafIntersection<GridImp>(center,facet)),
          leafFacet_(make_shared<MapType>())
    {
        if(facet==center->corners())
        {
            // This is the end facet.
            return;
        }

        for(std::size_t i=0; i<center->corners(); ++i)
            traverseAndPushLeafFacet(center->facet_[i], (*leafFacet_)[i]);


        // Search for the first intersection.
        // For each facet there is either one or it is a boundary intersection -> not anymore
        topLevelFacetIter_=leafFacet_->find(facet);
        assert(topLevelFacetIter_!=leafFacet_->end());
        assert(facet>0 || topLevelFacetIter_ == leafFacet_->begin());
        GridImp::getRealImplementation(intersection_).neighbor_=topLevelFacetIter_->second.begin();
        GridImp::getRealImplementation(intersection_).facetIndex_=facet;
        GridImp::getRealImplementation(intersection_).neighborEnd_=topLevelFacetIter_->second.end();
        assert(*GridImp::getRealImplementation(intersection_).neighbor_);

        if(center->facet_[facet]->elements_.size()==1){
            // This is a boundary facet.
            GridImp::getRealImplementation(intersection_).neighbor_= topLevelFacetIter_->second.end();
            return;
        }

        if(GridImp::getRealImplementation(intersection_).neighbor_ != GridImp::getRealImplementation(intersection_).neighborEnd_
           && *GridImp::getRealImplementation(intersection_).neighbor_==center)
        {
            ++GridImp::getRealImplementation(intersection_).neighbor_;
        }
    }

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLeafIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center)
        : intersection_(FoamGridLeafIntersection<GridImp>(center,center->corners())),
          leafFacet_(make_shared<MapType>())
    {
    }

public:

    typedef Dune::Intersection<const GridImp, typename Dune::FoamGridLeafIntersection<GridImp> > Intersection;

    //! equality
    bool equals(const FoamGridLeafIntersectionIterator<GridImp>& other) const {
        return GridImp::getRealImplementation(intersection_).equals(GridImp::getRealImplementation(other.intersection_));
    }

    //! prefix increment
    void increment()
    {
        if(GridImp::getRealImplementation(intersection_).facetIndex_ == GridImp::getRealImplementation(intersection_).center_->corners())
        {
            // This is already the end iterator
            DUNE_THROW(InvalidStateException, "Cannot increment a one past the end iterator");
            return;
        }
        assert(topLevelFacetIter_!=leafFacet_->end());

        if(GridImp::getRealImplementation(intersection_).neighbor_ !=
           GridImp::getRealImplementation(intersection_).neighborEnd_)
            // increment
            ++GridImp::getRealImplementation(intersection_).neighbor_;

        if(GridImp::getRealImplementation(intersection_).neighbor_ !=
           GridImp::getRealImplementation(intersection_).neighborEnd_ &&
           *(GridImp::getRealImplementation(intersection_).neighbor_) ==
           GridImp::getRealImplementation(intersection_).center_)
            // pointing to to the center_ -> increment another time
            ++GridImp::getRealImplementation(intersection_).neighbor_;

        if(GridImp::getRealImplementation(intersection_).neighbor_ ==
           GridImp::getRealImplementation(intersection_).neighborEnd_)
        {
            assert(topLevelFacetIter_!=leafFacet_->end());
            ++topLevelFacetIter_;
            ++GridImp::getRealImplementation(intersection_).facetIndex_;
            if(topLevelFacetIter_==leafFacet_->end())
                return;

            GridImp::getRealImplementation(intersection_).neighbor_
                = topLevelFacetIter_->second.begin();
            GridImp::getRealImplementation(intersection_).neighborEnd_=
                topLevelFacetIter_->second.end();
            if(GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.size()==1){
                // This is a boundary facet.
                GridImp::getRealImplementation(intersection_).neighbor_=
                    GridImp::getRealImplementation(intersection_).neighborEnd_;
                return;

            }
            if(GridImp::getRealImplementation(intersection_).neighbor_ != GridImp::getRealImplementation(intersection_).neighborEnd_
               && *GridImp::getRealImplementation(intersection_).neighbor_==
               GridImp::getRealImplementation(intersection_).center_)
            {
                ++GridImp::getRealImplementation(intersection_).neighbor_;
            }
        }
    }


    //! \brief dereferencing
    const Intersection & dereference() const {
        return intersection_;
    }

private:

    //! \brief Pushes all leaf facets into leafFacet_
    //!
    //! On returning leafFacet_ will contain the children
    //! of facet that do not have any children.
    //! \param facet The facet whose leafEdge we need.
    //!
    void traverseAndPushLeafFacet(FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* facet,
                                  ElementVector& leafFacet)
    {
        if(facet->isLeaf())
        {
            typedef typename ElementVector::const_iterator iter;
            for(iter eIt = facet->elements_.begin(); eIt!=facet->elements_.end(); ++eIt)
                leafFacet.push_back(*eIt);
            //leafFacet.insert(leafFacet.end(), facet->elements_.begin(), facet->elements_.end());
#ifndef NDEBUG
            for(iter k=leafFacet.begin(); k!=leafFacet.end(); ++k)
                assert(*k);
#endif
        }
        else
        {
            for (std::size_t i = 0; i < dimgrid; ++i)
                traverseAndPushLeafFacet(facet->sons_[i], leafFacet);
        }
    }

    Intersection intersection_;

    //! \brief pointer to map from facet index onto the intersections associated with the facet
    //!
    //! This has to be pointer to prevent copying it during copying the iterator.
    //! Otherwise the stored iterators to its entries would be invalidated.
    shared_ptr<MapType> leafFacet_;

    //! \brief Iterator pointing to the elements of the current facet.
    typename std::map<int, ElementVector>::const_iterator topLevelFacetIter_;
};




//! \todo Please doc me !
template<class GridImp>
class FoamGridLevelIntersectionIterator
{

    enum { dimgrid  = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };

    // Only the codim-0 entity is allowed to call the constructors
    friend class FoamGridEntity<0, dimgrid, GridImp>;

    template<typename, typename, typename>
    friend class Dune::IntersectionIterator;

    FoamGridLevelIntersectionIterator()
    {}

    //! \brief Constructor for a given grid entity and a given neighbor
    //! \param center Pointer to the element where the iterator was created.
    //! \param facet The index of the facet to start the investigation.
    FoamGridLevelIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center, std::size_t facet)
        : intersection_(FoamGridLevelIntersection<GridImp>(center,facet))
    {
        if(facet==center->corners())
        {
            // This is the end iterator
            return;
        }

        if(center->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.size()==1)
        {
            // This is a boundary facet.
            GridImp::getRealImplementation(intersection_).neighborEnd_=
                GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.end();
            GridImp::getRealImplementation(intersection_).neighbor_ = GridImp::getRealImplementation(intersection_).neighborEnd_;
            return;
        }

        // Search for the first intersection with a same level neighbor
        while(GridImp::getRealImplementation(intersection_).facetIndex_ != center->corners()) // not an  end iterator
        {
            GridImp::getRealImplementation(intersection_).neighbor_=GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.begin();
            GridImp::getRealImplementation(intersection_).neighborEnd_=GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.end();
            while(GridImp::getRealImplementation(intersection_).neighbor_!=GridImp::getRealImplementation(intersection_).neighborEnd_ &&
                  (GridImp::getRealImplementation(intersection_).center_==*GridImp::getRealImplementation(intersection_).neighbor_
                   ||GridImp::getRealImplementation(intersection_).center_->level()!=(*GridImp::getRealImplementation(intersection_).neighbor_)->level()))
            {
                ++GridImp::getRealImplementation(intersection_).neighbor_;
            }
            if(GridImp::getRealImplementation(intersection_).neighbor_==GridImp::getRealImplementation(intersection_).neighborEnd_)
            {
                if(GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.size()==1)
                {
                    // This is a boundary intersection.
                     GridImp::getRealImplementation(intersection_).neighbor_=GridImp::getRealImplementation(intersection_).neighborEnd_;
                    return;
                }else
                    // No valid intersection found on this facet, move to next one.
                    ++GridImp::getRealImplementation(intersection_).facetIndex_;
            }else
                // intersection with another element found!
                break;
        }
    }

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLevelIntersectionIterator(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center)
        : intersection_(FoamGridLevelIntersection<GridImp>(center,center->corners()))
    {
    }

public:

    typedef Dune::Intersection<const GridImp, typename Dune::FoamGridLevelIntersection<GridImp> > Intersection;

  //! equality
  bool equals(const FoamGridLevelIntersectionIterator<GridImp>& other) const {
      return GridImp::getRealImplementation(this->intersection_).equals(GridImp::getRealImplementation(other.intersection_));
  }

    //! prefix increment
    void increment() {
        if(GridImp::getRealImplementation(intersection_).facetIndex_==
           GridImp::getRealImplementation(intersection_).center_->corners())
        {
            // This is already the end iterator
            return;
        }
        if(GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.size()==1)
        {
            // This was a boundary intersection.
            ++GridImp::getRealImplementation(intersection_).facetIndex_;
            if(GridImp::getRealImplementation(intersection_).facetIndex_ < GridImp::getRealImplementation(intersection_).center_->corners()){
                // There is another facet, initialize neighbor_ iterator.
                GridImp::getRealImplementation(intersection_).neighbor_=GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.begin();
                GridImp::getRealImplementation(intersection_).neighborEnd_=GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.end();;
            }
        }else
        {
            // Move to the next intersection of this facet
            ++GridImp::getRealImplementation(intersection_).neighbor_;
        }

        // Search for the first intersection with a same level neighbor
        while(GridImp::getRealImplementation(intersection_).facetIndex_ != GridImp::getRealImplementation(intersection_).center_->corners()) // still a valid facet
        {
            if(GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.size()==1)
            {
                // This is a boundary intersection.
                GridImp::getRealImplementation(intersection_).neighbor_=GridImp::getRealImplementation(intersection_).neighborEnd_;
                return;
            }

            while(GridImp::getRealImplementation(intersection_).neighbor_!=GridImp::getRealImplementation(intersection_).neighborEnd_ &&
                  (GridImp::getRealImplementation(intersection_).center_==*GridImp::getRealImplementation(intersection_).neighbor_
                   ||GridImp::getRealImplementation(intersection_).center_->level()!=(*GridImp::getRealImplementation(intersection_).neighbor_)->level()))
            {
                // Wrong level or neighbor points to center. In both cases this intersection is invalid.
                ++GridImp::getRealImplementation(intersection_).neighbor_;
            }
            if(GridImp::getRealImplementation(intersection_).neighbor_==
               GridImp::getRealImplementation(intersection_).neighborEnd_)
            {
                if(GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.size()==1)
                {
                    // This is a boundary intersection.
                    GridImp::getRealImplementation(intersection_).neighbor_=GridImp::getRealImplementation(intersection_).neighborEnd_;
                    return;
                }
                else
                {
                    // No valid intersection found on this facet, move to next facet.
                    ++GridImp::getRealImplementation(intersection_).facetIndex_;
                    if(GridImp::getRealImplementation(intersection_).facetIndex_ < GridImp::getRealImplementation(intersection_).center_->corners())
                    {
                        // There is another facet, initialize neighbor_ iterator.
                        GridImp::getRealImplementation(intersection_).neighbor_=GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.begin();
                        GridImp::getRealImplementation(intersection_).neighborEnd_=GridImp::getRealImplementation(intersection_).center_->facet_[GridImp::getRealImplementation(intersection_).facetIndex_]->elements_.end();
                    }
                }
            }else
                // intersection with another element found!
                break;
        }
    }


    //! \brief dereferencing
    const Intersection & dereference() const
    {
        return intersection_;
    }
private:

    /** \brief The actual intersection
    */
    Intersection intersection_;

};


}  // namespace Dune

#endif
