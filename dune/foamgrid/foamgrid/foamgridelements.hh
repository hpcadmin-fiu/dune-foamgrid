#ifndef DUNE_FOAMGRID_ELEMENTS_HH
#define DUNE_FOAMGRID_ELEMENTS_HH

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridedge.hh>
#include <dune/common/nullptr.hh>
namespace Dune {

    template <int dimworld>
    class FoamGridEntityImp<2,dimworld>
        : public FoamGridEntityBase
    {
    public:

        /** \brief The different ways to mark an element for grid changes */
        enum MarkState { DO_NOTHING , COARSEN , REFINE };

        FoamGridEntityImp(int level, unsigned int id) 
            : FoamGridEntityBase(level,id),
              nSons_(0),
              markState_(DO_NOTHING), isNew_(false)
        {
          sons_[0]= sons_[1] = sons_[2] = sons_[3] = nullptr;
        }

        bool isLeaf() const {
            // The sons are either all nullptr or all != nullptr
            return sons_[0] == nullptr;
        }

        unsigned int corners() const {
            return 3;
        }

        GeometryType type() const {
            return GeometryType(GeometryType::simplex, 2);
        }

        /** \todo Implement me! */
        unsigned int nSons() const {
            return nSons_;
        }
      
        unsigned int nSons_;      

        array<FoamGridEntityImp<2,dimworld>*, 4> sons_;

        FoamGridEntityImp<2,dimworld>* father_;

        array<FoamGridEntityImp<1,dimworld>*, 3> edges_;

        FoamGridEntityImp<0,dimworld>* vertex_[3];
        
        /** \brief Stores requests for refinement and coarsening */
        MarkState markState_;
        
        /** \brief This flag is set by adapt() if this element has been newly created. */
        bool isNew_;
        
    };

}

#endif
