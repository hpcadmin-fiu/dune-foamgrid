#ifndef DUNE_FOAMGRID_LEVELITERATOR_HH
#define DUNE_FOAMGRID_LEVELITERATOR_HH

/** \file
* \brief The FoamGridLevelIterator class
*/

namespace Dune {

//**********************************************************************
//
// --FoamGridLevelIterator
/** \brief Iterator over all entities of a given codimension and level of a grid.
* \ingroup FoamGrid
*/
template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLevelIterator
{
    enum {dimgrid  = GridImp::dimension};
    enum {dimworld = GridImp::dimensionworld};

    using EntityImp = FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>;

public:

    using Entity = typename GridImp::template Codim<codim>::Entity;
    enum { codimension = codim };

    //! Constructor
    explicit FoamGridLevelIterator(const typename std::list<EntityImp>::const_iterator& it)
    : levelIterator_(it)
    {
        virtualEntity_.impl().setToTarget(&(*levelIterator_));
    }

    //! prefix increment
    void increment() {
        ++levelIterator_;
        virtualEntity_.impl().setToTarget(&(*levelIterator_));
    }

    //! dereferencing
    const Entity& dereference() const { return virtualEntity_; }

    //! equality
    bool equals(const FoamGridLevelIterator<codim, pitype, GridImp>& other) const {
      return virtualEntity_ == other.virtualEntity_;
    }


private:
    //! The entity that the iterator is pointing to
    Entity virtualEntity_;

    // This iterator derives from FoamGridEntityPointer, and that base class stores the value
    // of the iterator, i.e. the 'pointer' to the entity.  However, that pointer can not be
    // set to its successor in the level std::list, not even by magic.  Therefore we keep the
    // same information redundantly in this iterator, which can be incremented.
    typename std::list<EntityImp>::const_iterator levelIterator_;

};


}  // namespace Dune

#endif
