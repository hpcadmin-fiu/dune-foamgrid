// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_ELEMENTS_HH
#define DUNE_FOAMGRID_ELEMENTS_HH

#include <dune/common/fmatrix.hh>

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridedge.hh>
#include <dune/common/nullptr.hh>

namespace Dune {

    /** \brief Element specialization of FoamGridEntityImp. Element is a grid entity of topological codimension 0 and dimension dimgrid.*/
    template <int dimgrid, int dimworld>
    class FoamGridEntityImp<dimgrid, dimgrid, dimworld> {};

    /** \brief Element specialization of FoamGridEntityImp for 1d grids. Element is a grid entity of topological codimension 0 and dimension dimgrid.*/
    template <int dimworld>
    class FoamGridEntityImp<1, 1, dimworld>
        : public FoamGridEntityBase
    {
        /** \brief Grid dimension */
        enum {dimgrid = 1};

    public:

        /** \brief The different ways to mark an element for grid changes */
        enum MarkState { DO_NOTHING , COARSEN , REFINE, IS_COARSENED, SPAWN, PRUNE };

        FoamGridEntityImp(FoamGridEntityImp<0, dimgrid, dimworld>* v0,
                          FoamGridEntityImp<0, dimgrid, dimworld>* v1,
                          int level, unsigned int id)
            : FoamGridEntityBase(level,id), nSons_(0), father_(nullptr),
              refinementIndex_(-1), markState_(DO_NOTHING), isNew_(false), coarseningBlocked_(false)
        {
            vertex_[0] = v0;
            vertex_[1] = v1;
            facet_[0] = v0;
            facet_[1] = v1;
            sons_[0] = sons_[1] = nullptr;
        }


        FoamGridEntityImp(FoamGridEntityImp<0, dimgrid, dimworld>* v0,
                          FoamGridEntityImp<0, dimgrid, dimworld>* v1,
                          int level, unsigned int id,
                          FoamGridEntityImp* father)

            : FoamGridEntityBase(level,id), nSons_(0), father_(father),
              refinementIndex_(-1), markState_(DO_NOTHING), isNew_(false), coarseningBlocked_(false)
        {
            vertex_[0] = v0;
            vertex_[1] = v1;
            facet_[0] = v0;
            facet_[1] = v1;
            sons_[0] =sons_[1] = nullptr;
        }

        FoamGridEntityImp(int level, unsigned int id)
            : FoamGridEntityBase(level, id),
              nSons_(0), refinementIndex_(-1),
              markState_(DO_NOTHING), isNew_(false), coarseningBlocked_(false)
        {
          sons_[0] = sons_[1] = nullptr;
          father_ = nullptr;
        }

        /** \todo Implement this method! */
        bool isLeaf() const {
            return sons_[0]==nullptr;
        }

        /** \todo Implement me! */
        unsigned int nSons() const {
            return nSons_;
        }

        bool mightVanish() const
        {
            return markState_==COARSEN;
        }

        bool isNew() const
        {
            return isNew_;
        }

        GeometryType type() const {
            return GeometryType(GeometryType::simplex, 1);
        }

        bool hasFather() const
        {
            return father_!=nullptr;
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

        /** \brief Return level index of sub entity with codim = cc and local number i
         */
        int subLevelIndex (int i, unsigned int codim) const {
            assert(0<=codim && codim<=1);
            switch (codim) {
            case 0:
                return this->levelIndex_;
            case 1:
                return vertex_[i]->levelIndex_;
            }
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        /** \brief Return leaf index of sub entity with codim = cc and local number i
         */
        int subLeafIndex (int i,unsigned int codim) const {
            assert(0<=codim && codim<=1);
            switch (codim) {
            case 0:
                return this->leafIndex_;
            case 1:
                return vertex_[i]->leafIndex_;
            }
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        int refinementIndex_;

        MarkState markState_;

        /** \brief Stores the spawn point coordinates */
        Dune::FieldVector<double, dimworld> spawnPoint_;

        FoamGridEntityImp<0, dimgrid, dimworld>* vertex_[2];

        array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, 2> facet_;

        /** \brief links to refinements of this edge */
        array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 2> sons_;

        /** \brief The number of refined edges (0 or 2). */
        unsigned int nSons_;

        /** \brief Pointer to father element */
        FoamGridEntityImp<dimgrid, dimgrid, dimworld>* father_;

        /** \brief This flag is set by adapt() if this element has been newly created. */
        bool isNew_;

        /** \brief This flag is set by postGrow() if the element looses its right to coarsen
        *         because it contains a bifurcation facet without father
        */
        bool coarseningBlocked_;
    };

    /** \brief Element specialization of FoamGridEntityImp for 2d grids. Element is a grid entity of topological codimension 0 and dimension dimgrid.*/
    template <int dimworld>
    class FoamGridEntityImp<2, 2, dimworld>
        : public FoamGridEntityBase
    {
        /** \brief Grid dimension */
        enum {dimgrid = 2};

    public:

         /** \brief The different ways to mark an element for grid changes */
        enum MarkState { DO_NOTHING , COARSEN , REFINE, IS_COARSENED, SPAWN, PRUNE };

        FoamGridEntityImp(int level, unsigned int id)
            : FoamGridEntityBase(level,id),
              refinementIndex_(-1),
              nSons_(0),
              markState_(DO_NOTHING), isNew_(false), coarseningBlocked_(false)
        {
          sons_[0]= sons_[1] = sons_[2] = sons_[3] = nullptr;
          father_ = nullptr;
        }


        unsigned int corners() const {
            return 3;
        }

        GeometryType type() const {
            return GeometryType(GeometryType::simplex, 2);
        }


        bool hasFather() const
        {
            return father_!=nullptr;
        }

        bool mightVanish() const
        {
            return markState_==COARSEN;
        }

        bool isLeaf() const {
            // The sons are either all nullptr or all != nullptr
            return sons_[0] == nullptr;
        }

        bool isNew() const
        {
            return isNew_;
        }

        /** \todo Implement me! */
        unsigned int nSons() const {
            return nSons_;
        }

        /** \brief Compute local cordinates from global ones.
         * \param coord The global coordinates.
         * \return The corresponding local coordinates within the element.
         */

        FieldVector<double,2> globalToLocal(const FieldVector<double, dimworld>& coord) const
        {
            // If we set up the overdetermined system matrix we have
            // A[i][0]=vertex_[1].pos_[i]-vertex_[0].pos_[i];
            // A[i][1]=vertex_[2].pos_[i]-vertex_[0].pos_[i];
            // t[i]=coord[i]-vertex_[0].pos_[i];
            //
            // to determine the local coordinates we solve
            // A'A x= A' t
            //

            FieldMatrix<double,2,2> mat; // A'A
            // mat_{ij}=\sum_k A_{ki}A_{kj}
            mat=0;
            for(std::size_t i=0; i <dimworld; ++i)
            {
                mat[0][0]+=(vertex_[1]->pos_[i]-vertex_[0]->pos_[i])*(vertex_[1]->pos_[i]-vertex_[0]->pos_[i]);
            }
            for(std::size_t i=0; i <dimworld; ++i)
            {
                mat[1][0]+=(vertex_[2]->pos_[i]-vertex_[0]->pos_[i])*(vertex_[1]->pos_[i]-vertex_[0]->pos_[i]);
            }
            mat[0][1]=mat[1][0];
            for(std::size_t i=0; i <dimworld; ++i)
            {
                mat[1][1]+=(vertex_[2]->pos_[i]-vertex_[0]->pos_[i])*(vertex_[2]->pos_[i]-vertex_[0]->pos_[i]);
            }

            FieldVector<double, 2> b, x;
            b=0;
            for(std::size_t i=0; i <dimworld; ++i)
            {
                b[0]+=(vertex_[1]->pos_[i]-vertex_[0]->pos_[i])*(coord[i]-vertex_[0]->pos_[i]);
                b[1]+=(vertex_[2]->pos_[i]-vertex_[0]->pos_[i])*(coord[i]-vertex_[0]->pos_[i]);
            }
            mat.solve(x, b);
#ifndef NDEBUG
            FieldVector<double, dimworld> test(vertex_[0]->pos_);
            test.axpy(x[0], vertex_[1]->pos_);
            test.axpy(-x[0], vertex_[0]->pos_);
            test.axpy(x[1], vertex_[2]->pos_);
            test.axpy(-x[1], vertex_[0]->pos_);
            assert((test-coord).two_norm()< std::numeric_limits<double>::epsilon()*8);
#endif
            return x;
        }

        /** \brief Return level index of sub entity with codim = cc and local number i
         */
        int subLevelIndex (int i, unsigned int codim) const {
            assert(0<=codim && codim<=2);
            switch (codim) {
            case 0:
                return this->levelIndex_;
            case 1:
                return facet_[i]->levelIndex_;
            case 2:
                return vertex_[i]->levelIndex_;
            }
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        /** \brief Return leaf index of sub entity with codim = cc and local number i
         */
        int subLeafIndex (int i,unsigned int codim) const {
            assert(0<=codim && codim<=2);
            switch (codim) {
            case 0:
                return this->leafIndex_;
            case 1:
                return facet_[i]->leafIndex_;
            case 2:
                return vertex_[i]->leafIndex_;
            }
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        /**
         * \brief index of the refined element in the father
         *
         * For red refinement this is either the index of corner,
         * that is also a corner in the father element, within the father
         * or 3 if no corner is also a corner in the father.
         */
        int refinementIndex_;

        unsigned int nSons_;

        array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 4> sons_;

        FoamGridEntityImp<dimgrid, dimgrid ,dimworld>* father_;

        array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, 3> facet_;

        FoamGridEntityImp<0, dimgrid, dimworld>* vertex_[3];

        /** \brief Stores requests for refinement and coarsening */
        MarkState markState_;

        /** \brief Stores the spawn point or coordinates of the point to be pruned */
        Dune::FieldVector<double, dimworld> spawnPoint_;

        /** \brief This flag is set by adapt() if this element has been newly created. */
        bool isNew_;

        /** \brief This flag is set by postGrow() if the element looses its right to coarsen
        *         because it contains a bifurcation facet without father
        */
        bool coarseningBlocked_;

    };
}

#endif