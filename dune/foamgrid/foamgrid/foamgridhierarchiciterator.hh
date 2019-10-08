#ifndef DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH
#define DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH

/** \file
* \brief The FoamGridHierarchicIterator class
*/

#include <stack>

namespace Dune {


//**********************************************************************
//
/** \brief Iterator over the descendants of an entity.
* \ingroup FoamGrid
Mesh entities of codimension 0 ("elements") allow to visit all entities of
codimension 0 obtained through nested, hierarchic refinement of the entity.
Iteration over this set of entities is provided by the HierarchicIterator,
starting from a given entity.
*/
template<class GridImp>
class FoamGridHierarchicIterator
{
    enum {dimworld = GridImp::dimensionworld};
    enum {dimgrid  = GridImp::dimension};

    friend class FoamGridEntity<0, dimgrid,GridImp>;

    using StackEntry = const FoamGridEntityImp<dimgrid, dimgrid, dimworld, typename GridImp::ctype>*;

public:
    using Entity = typename GridImp::template Codim<0>::Entity;

    //! We only iterate over elements with this iterator
    enum { codimension = 0 };

    //! Constructor
    FoamGridHierarchicIterator(int maxlevel)
    : maxlevel_(maxlevel)
    {
        virtualEntity_.impl().setToTarget(nullptr);
    }

    //! \todo Please doc me !
    void increment()
    {
        if (elemStack.empty())
            return;

        StackEntry old_target = elemStack.top();
        elemStack.pop();

        // Traverse the tree no deeper than maxlevel
        if (old_target->level_ < maxlevel_) {

            // Load sons of old target onto the iterator stack
            if (!old_target->isLeaf())
                for (size_t i=0; i<old_target->nSons(); i++)
                    elemStack.push(old_target->sons_[i]);
        }

        virtualEntity_.impl().setToTarget(elemStack.empty() ?
                                          nullptr : elemStack.top());
    }

    //! dereferencing
    const Entity& dereference() const { return virtualEntity_; }

    //! equality
    bool equals(const FoamGridHierarchicIterator<GridImp>& other) const {
      return virtualEntity_ == other.virtualEntity_;
    }

private:
    //! The entity that the iterator is pointing to
    Entity virtualEntity_;

    //! max level to go down
    int maxlevel_;

    /** \brief For depth-first search */
    std::stack<StackEntry> elemStack;
};


}  // end namespace Dune

#endif
