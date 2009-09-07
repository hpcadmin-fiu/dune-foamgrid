#ifndef DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH
#define DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH

/** \file
* \brief The FoamGridHierarchicIterator class
*/

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
class FoamGridHierarchicIterator :
        public Dune::FoamGridEntityPointer <0,GridImp>
{
    public:
        
        typedef typename GridImp::template Codim<0>::Entity Entity;

#if 0    
        //! the default Constructor
        explicit FoamGridHierarchicIterator(const GridImp* identityGrid, const FoamGridElement& startEntity, int maxLevel) :
            FoamGridEntityPointer<0,GridImp>(identityGrid, startEntity.hostEntity_->hbegin(maxLevel)),
            identityGrid_(identityGrid),
            hostGridHierarchicIterator_(startEntity.hostEntity_->hbegin(maxLevel)),
            hostGridHierarchicEndIterator_(startEntity.hostEntity_->hend(maxLevel))
        {
            this->virtualEntity_.setToTarget(hostGridHierarchicIterator_);
        }
        
        
        //! \todo Please doc me !
        explicit FoamGridHierarchicIterator(const GridImp* identityGrid, const FoamGridElement& startEntity, int maxLevel, bool endDummy) :
            FoamGridEntityPointer<0,GridImp>(identityGrid, startEntity.hostEntity_->hend(maxLevel)),
            identityGrid_(identityGrid),
            hostGridHierarchicIterator_(startEntity.hostEntity_->hbegin(maxLevel)),
            hostGridHierarchicEndIterator_(startEntity.hostEntity_->hend(maxLevel))
    {}
    #endif
        
        //! \todo Please doc me !
        void increment()
        {
            if (elemStack.empty())
                return;
            
            const FoamGridElement* old_target = elemStack.top();
            elemStack.pop();
            
            // Traverse the tree no deeper than maxlevel
            if (old_target->level_ < maxlevel_) {
                
                // Load sons of old target onto the iterator stack
                if (!old_target->isLeaf()) {
                    
                    for (int i=0; i<old_target->nSons(); i++)
                        elemStack.push(old_target->sons_[i]);
                    
                }
                
            }
            
            this->virtualEntity_.setToTarget((elemStack.empty()) 
                                             ? NULL : elemStack.top());
        }

        
private:
        
    //! max level to go down 
    int maxlevel_;

    /** \brief For depth-first search */
    std::stack<const FoamGridElement*> elemStack;
};


}  // end namespace Dune

#endif
