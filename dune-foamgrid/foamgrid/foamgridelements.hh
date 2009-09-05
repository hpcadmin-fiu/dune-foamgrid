#ifndef DUNE_FOAMGRID_ELEMENTS_HH
#define DUNE_FOAMGRID_ELEMENTS_HH

#include <dune-foamgrid/foamgrid/foamgridvertex.hh>
#include <dune-foamgrid/foamgrid/foamgridedge.hh>

namespace Dune {

    template <int dimworld>
    class FoamGridEntityImp<2,dimworld>
    {
    public:

        /** \brief The different ways to mark an element for grid changes */
        enum MarkState { DO_NOTHING , COARSEN , REFINE };

        FoamGridEntityImp(int level, unsigned int id) 
            : id_(id), level_(level), 
              markState_(DO_NOTHING), isNew_(false)
        {
            sons_[0] = sons_[1] = NULL;
        }

        bool isLeaf() const {
            // The sons are either all NULL or all != NULL
            return sons_[0] == NULL;
        }

        unsigned int corners() const {
            return 3;
        }

        GeometryType type() const {
            return GeometryType(GeometryType::simplex, 2);
        }

        array<FoamGridEntityImp<2,dimworld>*, 4> sons_;

        FoamGridEntityImp<2,dimworld>* father_;

        array<FoamGridEntityImp<1,dimworld>*, 3> edges_;

        FoamGridEntityImp<0,dimworld>* vertex_[3];
        
        //! element number 
        unsigned int levelIndex_;
        
        unsigned int leafIndex_;
        
        /** \brief Unique and persistent id for elements */
        unsigned int id_;
        
        //! the level of the entity
        int level_;
        
        /** \brief Stores requests for refinement and coarsening */
        MarkState markState_;
        
        /** \brief This flag is set by adapt() if this element has been newly created. */
        bool isNew_;
        
    };

    typedef FoamGridEntityImp<2,3> FoamGridElement;

}

#endif
