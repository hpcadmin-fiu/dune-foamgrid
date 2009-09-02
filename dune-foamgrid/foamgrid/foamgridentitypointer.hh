#ifndef DUNE_IDENTITYGRID_ENTITY_POINTER_HH
#define DUNE_IDENTITYGRID_ENTITY_POINTER_HH

/** \file
* \brief The FoamGridEntityPointer class
*/

namespace Dune {


/** Acts as a pointer to an  entities of a given codimension.
*/
template<int codim, class GridImp>
class FoamGridEntityPointer
{
    private:
    
    enum { dim = GridImp::dimension };
    
    typedef typename SelectType<codim==0, FoamGridElement, FoamGridVertex>::Type TargetType;

    
    public:
    
    //! export the type of the EntityPointer Implementation.
    //! Necessary for the typeconversion between Iterators and EntityPointer
    typedef FoamGridEntityPointer EntityPointerImp;

    /** \brief Codimension of entity pointed to */
    enum { codimension = codim };

        typedef typename GridImp::template Codim<codim>::Entity Entity;
        
        typedef FoamGridEntityPointer<codim,GridImp> Base;
    
    //! Constructor from a FoamGrid entity
    FoamGridEntityPointer (const FoamGridEntity<codim,dim,GridImp>& entity)
        : virtualEntity_(entity.target_) 
    {}

    FoamGridEntityPointer (const typename std::list<TargetType>::const_iterator& it)
        : virtualEntity_(&(*it))
    {}

    FoamGridEntityPointer (const TargetType* target)
        : virtualEntity_(target)
    {}
        
        //! equality
        bool equals(const FoamGridEntityPointer<codim,GridImp>& i) const {
            return virtualEntity_.getTarget() == i.virtualEntity_.getTarget();
        }
    
        
        //! dereferencing
        Entity& dereference() const {
            return virtualEntity_;
        }

    //! Make this pointer as small as possible
    void compactify () {
        //virtualEntity_.getTarget().compactify();
    }
        
        //! ask for level of entity
        int level () const {
            return virtualEntity_.level();
        }
    
        
    protected:
    
        //! virtual entity
        mutable FoamGridMakeableEntity<codim,dim,GridImp> virtualEntity_;
    
        
};


} // end namespace Dune

#endif
