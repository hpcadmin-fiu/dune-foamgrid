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
    private:
    
        enum {dim = GridImp::dimension};
    
        
    public:
    
        //! \todo Please doc me !
        explicit FoamGridLeafIterator(const GridImp* identityGrid) :
            FoamGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template leafbegin<codim>()),
            hostGridLeafIterator_(identityGrid->hostgrid_->template leafbegin<codim>()),
            hostGridLeafEndIterator_(identityGrid->hostgrid_->template leafend<codim>())
        {
            this->virtualEntity_.setToTarget(hostGridLeafIterator_);
        }
    
        
        /** \brief Constructor which create the end iterator
        *  \param endDummy Here only to distinguish it from the other constructor
        */
        explicit FoamGridLeafIterator(const GridImp* identityGrid, bool endDummy) :
            FoamGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template leafend<codim>()),
            hostGridLeafIterator_(identityGrid->hostgrid_->template leafbegin<codim>()),
            hostGridLeafEndIterator_(identityGrid->hostgrid_->template leafend<codim>())
        {
        }
        
    
        //! prefix increment
        void increment() {
            ++hostGridLeafIterator_;
            this->virtualEntity_.setToTarget(hostGridLeafIterator_);
        }
    
    
    private:
    
        // /////////////////////////////////////
        //   Data members
        // /////////////////////////////////////
    
        // LevelIterator to the equivalent entity in the host grid
        typedef typename GridImp::HostGridType::template Codim<codim>::LeafIterator HostGridLeafIterator;
        
        //! \todo Please doc me !
        HostGridLeafIterator hostGridLeafIterator_;
        
        //! \todo Please doc me !
        HostGridLeafIterator hostGridLeafEndIterator_;
        
};


}  // namespace Dune
  
#endif
