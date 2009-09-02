#ifndef DUNE_FOAMGRID_VERTEX_HH
#define DUNE_FOAMGRID_VERTEX_HH

namespace Dune {

    template <int dim, int dimworld>
    class FoamGridEntityImp
    {};

    template <int dimworld>
    class FoamGridEntityImp<0,dimworld>
    {
    public:
        
        FoamGridEntityImp(int level, const Dune::FieldVector<double,dimworld>& pos) 
            : pos_(pos), level_(level), 
              son_(NULL)
        {}
        
        FoamGridEntityImp(int level, const FieldVector<double, dimworld>& pos, unsigned int id) 
            : pos_(pos), id_(id), level_(level), 
              son_(NULL) 
        {}
        
        //private: 
        bool isLeaf() const {
            return son_==NULL;
        }
        
        FieldVector<double, dimworld> pos_;
        
        //! entity number 
        unsigned int levelIndex_;
        
        unsigned int leafIndex_;
        
        unsigned int id_;
        
        //! level
        int level_;
        
        //! Son vertex on the next finer grid
        FoamGridEntityImp<0,dimworld>* son_;
        
    };

    /** \todo Only here for a transition */
    typedef FoamGridEntityImp<0,3> FoamGridVertex;

}

#endif
