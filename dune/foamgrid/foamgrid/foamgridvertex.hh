// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_VERTEX_HH
#define DUNE_FOAMGRID_VERTEX_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>


namespace Dune {

    /** \brief Base class for FoamGrid entity implementation classes */
    class FoamGridEntityBase
    {
    public:
        FoamGridEntityBase(int level, unsigned int id)
            : level_(level), id_(id), willVanish_()
        {}

        unsigned int level() const {
            return level_;
        }

        //! level
        int level_;

        //! entity number
        unsigned int levelIndex_;

        unsigned int leafIndex_;

        unsigned int id_;
        //! \brief Whether this entity will vanish due to coarsening.
        bool willVanish_;
    };

    /**
     * \brief The actual entity implementation
     *
     * \tparam dimentity The dimension of this entity
     * \tparam dimgrid The dimension of the grid
     * \tparam dimworld The world diemnsion
     */
    template <int dimentity, int dimgrid, int dimworld>
    class FoamGridEntityImp
    {};

    /** \brief Vertex specialization of FoamGridEntityImp */
    template <int dimgrid, int dimworld>
    class FoamGridEntityImp<0, dimgrid, dimworld>
        : public FoamGridEntityBase
    {
    public:

        FoamGridEntityImp(int level, const FieldVector<double, dimworld>& pos, unsigned int id)
            : FoamGridEntityBase(level, id),
              pos_(pos), son_(nullptr)
        {}

        //private:
        bool isLeaf() const {
            return son_==nullptr;
        }

        GeometryType type() const {
            return GeometryType(0);
        }

        /** \brief Number of corners (==1) */
        unsigned int corners() const {
            return 1;
        }

        FieldVector<double, dimworld> corner(int i) const {
            return pos_;
        }

        PartitionType partitionType() const {
            return InteriorEntity;
        }

        /** \brief Return level index of sub entity with codim = cc and local number i
         */
        int subLevelIndex (int i, unsigned int codim) const {
            assert(codim==dimgrid);
            return this->levelIndex_;
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        /** \brief Return leaf index of sub entity with codim = cc and local number i
         */
        int subLeafIndex (int i,unsigned int codim) const {
            assert(codim==dimgrid);
            return this->leafIndex_;
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        FieldVector<double, dimworld> pos_;

        //! Son vertex on the next finer grid
        FoamGridEntityImp<0, dimgrid, dimworld>* son_;

    };

}

#endif
