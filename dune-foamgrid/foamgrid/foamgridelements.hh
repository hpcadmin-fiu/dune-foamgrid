#ifndef DUNE_FOAMGRID_ELEMENTS_HH
#define DUNE_FOAMGRID_ELEMENTS_HH

namespace Dune {

    class FoamGridElement
    {
    public:

        /** \brief The different ways to mark an element for grid changes */
        enum MarkState { DO_NOTHING , COARSEN , REFINE };

        FoamGridElement(int level, unsigned int id) 
            : id_(id), level_(level), 
              markState_(DO_NOTHING), isNew_(false)
        {
            sons_[0] = sons_[1] = NULL;
        }

        bool isLeaf() const {
            DUNE_THROW(NotImplemented, "isLeaf()");
        }

        array<FoamGridElement*, 4> sons_;

        FoamGridElement* father_;

        FoamGridVertex* vertex_[3];
        
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



}

#endif
