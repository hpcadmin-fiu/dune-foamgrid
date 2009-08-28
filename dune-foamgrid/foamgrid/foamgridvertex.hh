#ifndef DUNE_FOAMGRID_VERTEX_HH
#define DUNE_FOAMGRID_VERTEX_HH

namespace Dune {

    class FoamGridVertex
    {
        enum {dimworld = 3};

    public:
        
        FoamGridVertex(int level, const Dune::FieldVector<double,dimworld>& pos) 
            : pos_(pos), level_(level), 
              son_(NULL)
        {}
        
        FoamGridVertex(int level, const FieldVector<double, dimworld>& pos, unsigned int id) 
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
        FoamGridVertex* son_;
        
    };



}

#endif
