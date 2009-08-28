#ifndef DUNE_IDENTITYGRID_INDEXSETS_HH
#define DUNE_IDENTITYGRID_INDEXSETS_HH

/** \file
* \brief The index and id sets for the FoamGrid class
*/

#include <vector>

namespace Dune {

    /** \todo Take the index types from the host grid */
    template<class GridImp>
    class FoamGridLevelIndexSet :
        public IndexSet<GridImp,FoamGridLevelIndexSet<GridImp> >
    {
    public:
        
        //! get index of an entity
        template<int codim>
        int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
        {
            return grid_->getRealImplementation(e).levelIndex();
        }
        
        
        //! get index of subEntity of a codim 0 entity
        int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const
        {
            return grid_->getRealImplementation(e).subLevelIndex(i,codim);
        }
        
        
        //! get number of entities of given codim, type and on this level
        int size (int codim) const {
            return grid_->size(level_,codim);
        }
        
        
        //! get number of entities of given codim, type and on this level
        int size (GeometryType type) const
        {
            return grid_->size(level_,type);
        }
        
        
        /** \brief Deliver all geometry types used in this grid */
        const std::vector<GeometryType>& geomTypes (int codim) const
        {
            return myTypes_[codim];
        }

        /** \brief Return true if the given entity is contained in the index set */
        template<class EntityType>
        bool contains (const EntityType& e) const
        {
            DUNE_THROW(NotImplemented, "contains()");
        }
        
        /** \brief Set up the index set */
        void update(const GridImp& grid, int level)
        {
            grid_  = &grid;
            level_ = level;

            // ///////////////////////////////
            //   Init the element indices
            // ///////////////////////////////
            numElements_ = 0;
            std::list<FoamGridElement>::const_iterator eIt;
            for (eIt = grid_->elements_[level_].begin(); eIt!=grid_->elements_[level_].end(); ++eIt)
             /** \todo Remove this const cast */
                 *const_cast<unsigned int*>(&(eIt->levelIndex_)) = numElements_++;
            
            // //////////////////////////////
            //   Init the vertex indices
            // //////////////////////////////
            
            numVertices_ = 0;
            std::list<FoamGridVertex>::const_iterator vIt;
            for (vIt = grid_->vertices_[level_].begin(); vIt!=grid_->vertices_[level_].end(); ++vIt)
                /** \todo Remove this const cast */
                *const_cast<unsigned int*>(&(vIt->levelIndex_)) = numVertices_++;
            
            // ///////////////////////////////////////////////
            //   Update the list of geometry types present
            // ///////////////////////////////////////////////
            if (numElements_>0) {
                myTypes_[0].resize(1);
                myTypes_[0][0] = GeometryType(1);
            } else
                myTypes_[0].resize(0);
            
            if (numVertices_>0) {
                myTypes_[1].resize(1);
                myTypes_[1][0] = GeometryType(0);
            } else
                myTypes_[1].resize(0);
        }
        
        
        GridImp* grid_;
        
        int level_;

        int numElements_;

        int numEdges_;

        int numVertices_;

        /** \brief The GeometryTypes present for each codim */
        std::vector<GeometryType> myTypes_[3];
    };


template<class GridImp>
class FoamGridLeafIndexSet :
    public IndexSet<GridImp,FoamGridLeafIndexSet<GridImp> >
{
    
public:
    
        //! constructor stores reference to a grid and level
        FoamGridLeafIndexSet (const GridImp& grid)
            : grid_(&grid)
        {}
        
        
        //! get index of an entity
        /*
            We use the RemoveConst to extract the Type from the mutable class,
            because the const class is not instantiated yet.
        */
        template<int codim>
        int index (const typename remove_const<GridImp>::type::template Codim<codim>::Entity& e) const
        {
            return grid_.getRealImplementation(e).leafIndex(); 
        }
        
        
        //! get index of subEntity of a codim 0 entity
        /*
            We use the RemoveConst to extract the Type from the mutable class,
            because the const class is not instantiated yet.
        */
    int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, 
                  int i, 
                  int codim) const
    {
        return grid_.getRealImplementation(e).subLeafIndex(i,codim);
    }
        
        
    //! get number of entities of given type
    int size (GeometryType type) const
    {
        if (type.isVertex()) {

            return numVertices_;

        } else if (type.isLine()) {

            return numEdges_;

        } else if (type.dim()==2)

            return numElements_;

        return 0;
    }
        
        
        //! get number of entities of given codim
        int size (int codim) const
        {
            if (codim==2)
                return numVertices_;
            
            if (codim==1)
                return numEdges_;
          
            if (codim==0)
                return numElements_;

            return 0;
        }
        
        
        /** \brief Deliver all geometry types used in this grid */
        const std::vector<GeometryType>& geomTypes (int codim) const
        {
            return myTypes_[codim];
        }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
        DUNE_THROW(NotImplemented, "contains");
    }


        
    /** Recompute the leaf numbering */
    void update(const GridImp& grid)
    {
#if 0
        // ///////////////////////////////
        //   Init the element indices
        // ///////////////////////////////
        numElements_ = 0;
        typename GridImp::Traits::template Codim<0>::LeafIterator eIt    = grid_.template leafbegin<0>();
        typename GridImp::Traits::template Codim<0>::LeafIterator eEndIt = grid_.template leafend<0>();
        
        for (; eIt!=eEndIt; ++eIt)
            grid_.getRealImplementation(*eIt).target_->leafIndex_ = numElements_++;

        // //////////////////////////////
        //   Init the vertex indices
        // //////////////////////////////
        
        numVertices_ = 0;
        
        for (int i=grid_.maxLevel(); i>=0; i--) {

            std::list<FoamGridVertex>::const_iterator vIt;
            for (vIt = grid_.vertices[i].begin(); vIt!=grid_.vertices[i].end(); ++vIt) {
                
                /** \todo Remove the const casts */
                if (vIt->isLeaf())
                    const_cast<OneDEntityImp<0>*>(vIt)->leafIndex_ = numVertices_++;
                else
                    const_cast<OneDEntityImp<0>*>(vIt)->leafIndex_ = vIt->son_->leafIndex_;

            }

        }

        // ///////////////////////////////////////////////
        //   Update the list of geometry types present
        // ///////////////////////////////////////////////
        if (numElements_>0) {
            myTypes_[0].resize(1);
            myTypes_[0][0] = GeometryType(1);
        } else
            myTypes_[0].resize(0);

        if (numVertices_>0) {
            myTypes_[1].resize(1);
            myTypes_[1][0] = GeometryType(0);
        } else
            myTypes_[1].resize(0);
#endif
        }
   
        
        GridImp* grid_;

        int numElements_;

        int numEdges_;

        int numVertices_;

        /** \brief The GeometryTypes present for each codim */
        std::vector<GeometryType> myTypes_[3];

};




template <class GridImp>
class FoamGridGlobalIdSet :
    public IdSet<GridImp,FoamGridGlobalIdSet<GridImp>, unsigned int>
{
            
    public:
        //! constructor stores reference to a grid
        FoamGridGlobalIdSet (const GridImp& g) : grid_(&g) {}
        
        //! define the type used for persistent indices
        typedef unsigned int IdType;
        
        
        //! get id of an entity
        /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instantiated yet.
        */
        template<int cd>
        IdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
        {
            return grid_.getRealImplementation(e).globalId();
        }
    
        
        //! get id of subEntity
        /*
            We use the remove_const to extract the Type from the mutable class,
            because the const class is not instantiated yet.
        */
    IdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
        {
            return grid_.getRealImplementation(e).subId(i,codim);
        }

        
        /** \todo Should be private */
        void update() {}

        
        const GridImp* grid_;
};




template<class GridImp>
class FoamGridLocalIdSet :
    public IdSet<GridImp,FoamGridLocalIdSet<GridImp>, unsigned int>
{
        
    public:
        //! define the type used for persistent local ids
        typedef unsigned int IdType;

        
        //! constructor stores reference to a grid
        FoamGridLocalIdSet (const GridImp& g) : grid_(&g) {}
    
        
        //! get id of an entity
        /*
            We use the remove_const to extract the Type from the mutable class,
            because the const class is not instantiated yet.
        */
        template<int cd>
        IdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
        {
            // Return id of the host entity
            return grid_.getRealImplementation(e).globalId();
        }
        
        
        //! get id of subEntity
        /*
        * We use the remove_const to extract the Type from the mutable class,
        * because the const class is not instantiated yet.
        */
    IdType subId (const typename remove_const<GridImp>::type::template Codim<0>::Entity& e, int i, int codim) const
        {
            return grid_.getRealImplementation(e).subId(i,codim);
        }
        
        
        /** \todo Should be private */
        void update() {}

        
        const GridImp* grid_;
};


}  // namespace Dune


#endif
