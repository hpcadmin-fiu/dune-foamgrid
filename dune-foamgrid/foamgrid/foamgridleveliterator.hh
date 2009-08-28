#ifndef DUNE_FOAMGRID_LEVELITERATOR_HH
#define DUNE_FOAMGRID_LEVELITERATOR_HH

/** \file
* \brief The FoamGridLevelIterator class and its specializations
*/

namespace Dune {




//**********************************************************************
//
// --FoamGridLevelIterator
/** \brief Iterator over all entities of a given codimension and level of a grid.
* \ingroup FoamGrid
*/
template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLevelIterator :
    public Dune::FoamGridEntityPointer <codim,GridImp>
{
    private:

    typedef typename SelectType<codim==0, FoamGridElement, FoamGridVertex>::Type TargetType;

    public:
        
        //! Constructor
        explicit FoamGridLevelIterator(const typename std::list<TargetType>::const_iterator& it)
            : FoamGridEntityPointer<codim,GridImp>(it),
              levelIterator_(it)
        {
            this->virtualEntity_.setToTarget(&(*levelIterator_));
        }
        
    //! prefix increment
        void increment() {
            ++levelIterator_;
            this->virtualEntity_.setToTarget(&(*levelIterator_));
        }
        
        
    private:
    
        //! \todo Please doc me !
    typename std::list<TargetType>::const_iterator levelIterator_;
        
};


}  // namespace Dune
  
#endif
