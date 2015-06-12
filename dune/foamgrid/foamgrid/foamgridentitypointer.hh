#ifndef DUNE_FOAMGRID_ENTITY_POINTER_HH
#define DUNE_FOAMGRID_ENTITY_POINTER_HH

/** \file
* \brief The FoamGridEntityPointer class
*/

#include <list>

#include <dune/foamgrid/foamgrid/foamgridentity.hh>

namespace Dune {


/** Acts as a pointer to an  entities of a given codimension.
*/
template<int codim, class GridImp>
class FoamGridEntityPointer
{
    private:

    enum { dimgrid  = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };

    public:

    //! export the type of the EntityPointer Implementation.
    //! Necessary for the typeconversion between Iterators and EntityPointer
    typedef FoamGridEntityPointer EntityPointerImp;

    /** \brief Codimension of entity pointed to */
    enum { codimension = codim };

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    FoamGridEntityPointer()
      : virtualEntity_(FoamGridEntity<codim, dimgrid, GridImp>())
    {}

    //! Constructor from a FoamGrid entity
    FoamGridEntityPointer (const FoamGridEntity<codim, dimgrid, GridImp>& entity)
        : virtualEntity_(entity)
    {}

    FoamGridEntityPointer (const typename std::list<FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld> >::const_iterator& it)
        : virtualEntity_(FoamGridEntity<codim, dimgrid, GridImp>())
    {
        GridImp::getRealImplementation(virtualEntity_).setToTarget(&*it);
    }

    FoamGridEntityPointer (const FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>* it)
        : virtualEntity_(FoamGridEntity<codim, dimgrid, GridImp>())
    {
        GridImp::getRealImplementation(virtualEntity_).setToTarget(it);
    }
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,4)
        //! equality
        bool equals(const FoamGridEntityPointer<codim,GridImp>& other) const {
            return virtualEntity_ == other.virtualEntity_;
        }
#else
        bool equals(const FoamGridEntityPointer<codim,GridImp>& i) const {
            return GridImp::getRealImplementation(virtualEntity_).target_ == GridImp::getRealImplementation(i.virtualEntity_).target_;
        }
#endif

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,4)
        //! dereferencing
        const Entity& dereference() const {
            return virtualEntity_;
        }
#else
        //! dereferencing
        Entity& dereference() const {
            return virtualEntity_;
        }
#endif
        //! ask for level of entity
        int level () const {
            return virtualEntity_.level();
        }

    protected:
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,4)
        //! virtual entity
        Entity virtualEntity_;
#else
        mutable MakeableInterfaceObject<Entity> virtualEntity_;
#endif
};


} // end namespace Dune

#endif
