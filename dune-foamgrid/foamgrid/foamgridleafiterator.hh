#ifndef DUNE_IDENTITYGRIDLEAFITERATOR_HH
#define DUNE_IDENTITYGRIDLEAFITERATOR_HH

/** \file
* \brief The FoamGridLeafIterator class
*/

namespace Dune {


/** \brief Iterator over all entities of a given codimension and level of a grid.
*  \ingroup FoamGrid
*/
template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLeafIterator :
    public Dune::FoamGridEntityPointer <codim,GridImp>
{
    enum {dim = GridImp::dimension};

    //friend class OneDGridEntity<codim,dim,GridImp>;

    typedef typename SelectType<codim==0, FoamGridElement, FoamGridVertex>::Type TargetType;

public:

    FoamGridLeafIterator(const GridImp& grid) 
        : grid_(&grid),
          FoamGridEntityPointer <codim,GridImp>(NULL)
    {

        /** \todo Can a make the fullRefineLevel work somehow? */
        const int fullRefineLevel = 0;
        
        if (codim==0)
            // The &* turns an iterator into a plain pointer
            this->virtualEntity_.setToTarget((TargetType*)&*grid_->elements_[fullRefineLevel].begin());
        else
            this->virtualEntity_.setToTarget((TargetType*)&*grid_->vertices_[fullRefineLevel].begin());

        if (!this->virtualEntity_.getTarget()->isLeaf())
            increment();
    }

  //! Constructor
    FoamGridLeafIterator() 
        : FoamGridEntityPointer <codim,GridImp>(NULL),
          grid_(NULL),
          levelIterator_(grid_->elements_[0].begin())
    {}

    //! prefix increment
    void increment() {
        // Increment until you find a leaf entity
        do {
            globalIncrement();
        } while (this->virtualEntity_.getTarget() && !this->virtualEntity_.getTarget()->isLeaf());
    }

private:

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {

        // Backup current level because it may not be accessible anymore after
        // setting the pointer to the next entity.
        const int oldLevel = this->virtualEntity_.level();

        // Increment on this level
        ++levelIterator_;
        this->virtualEntity_.setToTarget(&(*levelIterator_));

        // If beyond the end of this level set to first of next level
        if (!this->virtualEntity_.getTarget() && oldLevel < grid_->maxLevel()) {

            if (codim==0) {
                // cast is necessary to make the code compile.  If this branch is taken the
                // cast is empty
                levelIterator_ = *(typename std::list<TargetType>::const_iterator*)&grid_->elements_[oldLevel+1].begin();
                this->virtualEntity_.setToTarget((TargetType*)&*grid_->elements_[oldLevel+1].begin());
            } else {
                // cast is necessary to make the code compile.  If this branch is taken the
                // cast is empty
                levelIterator_ = *(typename std::list<TargetType>::const_iterator*)&grid_->vertices_[oldLevel+1].begin();
                this->virtualEntity_.setToTarget((TargetType*)&*grid_->vertices_[oldLevel+1].begin());
            }

        }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

    //! \todo Please doc me !
    typename std::list<TargetType>::const_iterator levelIterator_;
};


}  // namespace Dune
  
#endif
