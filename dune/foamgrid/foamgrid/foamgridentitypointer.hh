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

    //! Constructor from a FoamGrid entity
    FoamGridEntityPointer (const FoamGridEntity<codim, dimgrid, GridImp>& entity)
        : virtualEntity_(entity.target_)
    {}

    FoamGridEntityPointer (const typename std::list<FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld> >::const_iterator& it)
        : virtualEntity_(&(*it))
    {}

    FoamGridEntityPointer (const FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>* target)
        : virtualEntity_(target)
    {}

        //! equality
        bool equals(const FoamGridEntityPointer<codim,GridImp>& other) const {
            return virtualEntity_ == other.virtualEntity_;
        }


        //! dereferencing
        Entity& dereference() const {
            return virtualEntity_;
        }

        //! ask for level of entity
        int level () const {
            return virtualEntity_.level();
        }

    protected:

        //! virtual entity
        Entity virtualEntity_;
};


} // end namespace Dune

#endif
