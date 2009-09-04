#ifndef DUNE_FOAMGRID_EDGE_HH
#define DUNE_FOAMGRID_EDGE_HH

#include <dune-foamgrid/foamgrid/foamgridvertex.hh>

namespace Dune {

    template <int dimworld>
    class FoamGridEntityImp<1,dimworld>
    {
    public:

        FoamGridEntityImp(int level, unsigned int id) 
            : id_(id), level_(level)
        {}

        bool isLeaf() const {
            DUNE_THROW(NotImplemented, "isLeaf()");
        }

        unsigned int boundaryId() const {
            return boundaryId_;
        }

        std::vector<FoamGridEntityImp<2,dimworld>*> elements_;

        FoamGridEntityImp<0,dimworld>* vertex_[2];

        //! element number 
        unsigned int levelIndex_;
        
        unsigned int leafIndex_;
        
        /** \brief Unique and persistent id for elements */
        unsigned int id_;

        /** \brief The boundary id.  Only used if this edge is a boundary edge */
        unsigned int boundaryId_;
        
        //! the level of the entity
        int level_;
        
    };

    typedef FoamGridEntityImp<1,3> FoamGridEdge;

}

#endif
