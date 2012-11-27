// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:

#ifndef DUNE_FOAMGRID_EDGE_HH
#define DUNE_FOAMGRID_EDGE_HH

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>

namespace Dune {

    template <int dimworld>
    class FoamGridEntityImp<1,dimworld>
        : public FoamGridEntityBase
    {
    public:

        FoamGridEntityImp(const FoamGridEntityImp<0,dimworld>* v0, 
                          const FoamGridEntityImp<0,dimworld>* v1, 
                          int level, unsigned int id) 
            : FoamGridEntityBase(level,id), nSons_(0)
        {
            vertex_[0] = v0;
            vertex_[1] = v1;
            sons_[0] =sons_[1] = nullptr;
        }

        /** \todo Implement this method! */
        bool isLeaf() const {
            return nSons_<2;
        }

        unsigned int boundarySegmentIndex() const {
            return boundaryId_;
        }

        GeometryType type() const {
            return GeometryType(1);
        }

        /** \brief Number of corners (==2) */
        unsigned int corners() const {
            return 2;
        }

        FieldVector<double, dimworld> corner(int i) const {
            return vertex_[i]->pos_;
        }

        PartitionType partitionType() const {
            return InteriorEntity;
        }

        std::vector<const FoamGridEntityImp<2,dimworld>*> elements_;

        const FoamGridEntityImp<0,dimworld>* vertex_[2];

        /** \brief The boundary id.  Only used if this edge is a boundary edge */
        unsigned int boundaryId_;

        /** \brief links to refinements of this edge */
        array<FoamGridEntityImp<1,dimworld>*,2> sons_;

        /** \brief The number of refined edges (0 or 2). */
        unsigned int nSons_;
        
    };

}

#endif
