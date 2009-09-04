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
class FoamGridLevelIterator :
    public Dune::FoamGridEntityPointer <codim,GridImp>
{
    enum {dim      = GridImp::dimension};
    enum {dimworld = GridImp::dimensionworld};

    public:
        
        //! Constructor
    explicit FoamGridLevelIterator(const typename std::list<FoamGridEntityImp<dim-codim,dimworld> >::const_iterator& it)
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
    typename std::list<FoamGridEntityImp<dim-codim,dimworld> >::const_iterator levelIterator_;
        
};


}  // namespace Dune
  
#endif
