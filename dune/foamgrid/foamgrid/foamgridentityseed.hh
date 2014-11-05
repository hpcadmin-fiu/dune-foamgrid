#ifndef DUNE_FOAMGRID_ENTITY_SEED_HH
#define DUNE_FOAMGRID_ENTITY_SEED_HH

/**
 * \file
 * \brief The SubGridEntitySeed class
 */

#include <dune/common/nullptr.hh>

#include <dune/foamgrid/foamgrid/foamgridentity.hh>

namespace Dune {


/**
 * \brief The EntitySeed class provides the minmal information needed to restore an Entity using the grid.
 * \ingroup SubGrid
 *
 */
template<int codim, class GridImp>
class FoamGridEntitySeed
{
        template<int dimgrid, int dimworld>
        friend class FoamGrid;

    protected:

        enum {dim = GridImp::dimension};
        enum {dimworld = GridImp::dimensionworld};
        enum {mydim = dim-codim};

        // Entity type of the hostgrid
        typedef FoamGridEntityImp<mydim, dimworld> EntityImplType;

    public:

        enum {codimension = codim};

        /**
         * \brief Create EntitySeed from hostgrid Entity
         *
         * We call hostEntity.seed() directly in the constructor
         * of SubGridEntitySeed to allow for return value optimization.
         *
         * If would use SubGridEntitySeed(hostEntity.seed())
         * we would have one copy even with optimization enabled.
         */
        FoamGridEntitySeed() :
            entityImplPointer_(nullptr)
        {}

        /**
         * \brief Create EntitySeed from hostgrid Entity
         *
         * We call hostEntity.seed() directly in the constructor
         * of SubGridEntitySeed to allow for return value optimization.
         *
         * If would use SubGridEntitySeed(hostEntity.seed())
         * we would have one copy even with optimization enabled.
         */
        FoamGridEntitySeed(const EntityImplType* impl) :
            entityImplPointer_(impl)
        {}

        /** \brief check whether it is safe to create an Entity from this Seed */
        bool isValid() const
        {
          return entityImplPointer_;
        }

    protected:

        const EntityImplType* getImplementationPointer() const
        {
            return entityImplPointer_;
        }

    private:

        const FoamGridEntityImp<mydim, dimworld>* entityImplPointer_;
};

} // namespace Dune


#endif
