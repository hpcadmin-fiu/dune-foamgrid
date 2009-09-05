#ifndef DUNE_FOAMGRID_HH
#define DUNE_FOAMGRID_HH

/** \file
* \brief The FoamGrid class
*/

#include <list>
#include <map>

#include <dune/common/collectivecommunication.hh>
#include <dune/common/tuples.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// Implementation classes
#include "foamgrid/foamgridvertex.hh"
#include "foamgrid/foamgridedge.hh"
#include "foamgrid/foamgridelements.hh"

// The components of the FoamGrid interface
#include "foamgrid/foamgridgeometry.hh"
#include "foamgrid/foamgridentity.hh"
#include "foamgrid/foamgridentitypointer.hh"
#include "foamgrid/foamgridintersectioniterators.hh"
#include "foamgrid/foamgridleveliterator.hh"
#include "foamgrid/foamgridleafiterator.hh"
#include "foamgrid/foamgridhierarchiciterator.hh"
#include "foamgrid/foamgridindexsets.hh"

namespace Dune {

// Forward declaration
class FoamGrid;

template<int codim>                        
class FoamGridLevelIteratorFactory;


template<int dim, int dimworld>
struct FoamGridFamily
{
    typedef GridTraits<
        dim,   // dim
        dimworld,   // dimworld
        Dune::FoamGrid,
        FoamGridGeometry,
        FoamGridEntity,
        FoamGridEntityPointer,
        FoamGridLevelIterator,
        FoamGridLeafIntersection,
        FoamGridLevelIntersection,
        FoamGridLeafIntersectionIterator,
        FoamGridLevelIntersectionIterator,
        FoamGridHierarchicIterator,
        FoamGridLeafIterator,
        FoamGridLevelIndexSet< const FoamGrid >,
        FoamGridLeafIndexSet< const FoamGrid >,
        FoamGridGlobalIdSet< const FoamGrid >,
        unsigned int,   // global id type
        FoamGridLocalIdSet< const FoamGrid >,
        unsigned int,   // local id type
        CollectiveCommunication<FoamGrid>
            > Traits;
};




//**********************************************************************
//
// --FoamGrid
//
//**********************************************************************

/** \brief [<em> provides \ref Dune::Grid </em>]
*
*/
class FoamGrid :
        public GridDefaultImplementation  <2, 3, double, FoamGridFamily<2,3> >
{
    
    friend class FoamGridLevelIteratorFactory <0>;
    friend class FoamGridLevelIteratorFactory <2>;

    friend class FoamGridLevelIndexSet<const FoamGrid >;
    friend class FoamGridLeafIndexSet<const FoamGrid >;
    friend class FoamGridGlobalIdSet<const FoamGrid >;
    friend class FoamGridLocalIdSet<const FoamGrid >;
    friend class FoamGridHierarchicIterator<const FoamGrid >;
    friend class FoamGridLevelIntersectionIterator<const FoamGrid >;
    friend class FoamGridLeafIntersectionIterator<const FoamGrid >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLeafIterator;
    
    template <class GridType_>
    friend class GridFactory;

    template<int codim_, int dim_, class GridImp_>
    friend class FoamGridEntity;

    public:
        
    //**********************************************************
    // The Interface Methods
    //**********************************************************
    
    //! type of the used GridFamily for this grid
    typedef FoamGridFamily<2,3>  GridFamily;
    
    //! the Traits
    typedef FoamGridFamily<2,3>::Traits Traits;
    
    //! The type used to store coordinates, inherited from the HostGrid
    typedef double ctype;
    
    /** \brief Constructor
     */
    FoamGrid() 
        : leafIndexSet_(*this),
          globalIdSet_(*this),
          localIdSet_(*this)
    {}
        
        //! Desctructor
        ~FoamGrid()
        {
            // Delete level index sets
            for (size_t i=0; i<levelIndexSets_.size(); i++)
                if (levelIndexSets_[i])
                    delete (levelIndexSets_[i]);
        }
        
        
        //! return grid name
        std::string name() const
        {
            return "FoamGrid";
        }
    
        
        //! Return maximum level defined in this grid. Levels are numbered
        //! 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const {
            return entityImps_.size()-1;;
        }
        
        
        //! Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid>(Dune::get<dimension-codim>(entityImps_[level]).begin());
        }
    
        
        //! one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");
            
            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid>(Dune::get<dimension-codim>(entityImps_[level]).end());
        }
        
        
        //! Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");
            
            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid>(Dune::get<dimension-codim>(entityImps_[level]).begin());
        }
        

        //! one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");
            
            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid>(Dune::get<dimension-codim>(entityImps_[level]).end());
        }
        
    
        //! Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
            return FoamGridLeafIterator<codim,All_Partition, const FoamGrid >(*this);
        }
        
    
        //! one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const {
            return FoamGridLeafIterator<codim,All_Partition, const FoamGrid >();
        }
        
    
        //! Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
            return FoamGridLeafIterator<codim,PiType, const FoamGrid >(*this);
        }
        
        
        //! one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
            return FoamGridLeafIterator<codim,PiType, const FoamGrid >();
        }
        

        /** \brief Number of grid entities per level and codim
         */
        int size (int level, int codim) const {

            // Turn dynamic index into static index
            if (codim==0)
                return Dune::get<0>(entityImps_[level]).size();
            if (codim==1)
                return Dune::get<1>(entityImps_[level]).size();
            if (codim==2)
                return Dune::get<2>(entityImps_[level]).size();

            return 0;
        }
        
        
        //! number of leaf entities per codim in this process
        int size (int codim) const{
            return leafIndexSet().size(codim);
        }
        
        
        //! number of entities per level, codim and geometry type in this process
        int size (int level, GeometryType type) const {
            return levelIndexSets_[level]->size(type);
        }
        
            
        //! number of leaf entities per codim and geometry type in this process
        int size (GeometryType type) const
        {
            return leafIndexSet().size(type);
        }
        
        
        /** \brief Access to the GlobalIdSet */
        const Traits::GlobalIdSet& globalIdSet() const{
            return globalIdSet_;
        }
        
        
        /** \brief Access to the LocalIdSet */
        const Traits::LocalIdSet& localIdSet() const{
            return localIdSet_;
        }
        
        
        /** \brief Access to the LevelIndexSets */
        const Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return *levelIndexSets_[level];
        }
        
        
        /** \brief Access to the LeafIndexSet */
        const Traits::LeafIndexSet& leafIndexSet() const
        {
            return leafIndexSet_;
        }
        
        
        /** @name Grid Refinement Methods */
        /*@{*/
        
        
        /** global refinement
        */
        void globalRefine (int refCount)
        {
            DUNE_THROW(NotImplemented, "globalRefine!");
        }
        
        /** \brief Mark entity for refinement
        *
        * This only works for entities of codim 0.
        * The parameter is currently ignored
        *
        * \return <ul>
        * <li> true, if marking was succesfull </li>
        * <li> false, if marking was not possible </li>
        * </ul>
        */
        bool mark(int refCount, const Traits::Codim<0>::EntityPointer & e)
        {
            return false;
        }
        
        /** \brief Return refinement mark for entity
        *
        * \return refinement mark (1,0,-1)
        */
        int getMark(const Traits::Codim<0>::EntityPointer & e) const
        {
            return 0;
        }

        //! \todo Please doc me !
        bool preAdapt() {
            DUNE_THROW(NotImplemented, "preAdapt");
        }
        
        
        //! Triggers the grid refinement process
        bool adapt()
        {
            DUNE_THROW(NotImplemented, "adapt");
        }

        /** \brief Clean up refinement markers */
        void postAdapt() {
            DUNE_THROW(NotImplemented, "postAdapt");
        }
        
        /*@}*/
        
        /** \brief Size of the overlap on the leaf level */
        unsigned int overlapSize(int codim) const {
            return 0;
        }
        
        
        /** \brief Size of the ghost cell layer on the leaf level */
        unsigned int ghostSize(int codim) const {
            return 0;
        }
        
        
        /** \brief Size of the overlap on a given level */
        unsigned int overlapSize(int level, int codim) const {
            return 0;
        }
        
        
        /** \brief Size of the ghost cell layer on a given level */
        unsigned int ghostSize(int level, int codim) const {
            return 0;
        }
        
            
#if 0
        /** \brief Distributes this grid over the available nodes in a distributed machine
        *
        * \param minlevel The coarsest grid level that gets distributed
        * \param maxlevel does currently get ignored
        */
        void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
            DUNE_THROW(NotImplemented, "FoamGrid::loadBalance()");
        }
        
        /** \brief The communication interface
        *  @param T: array class holding data associated with the entities
        *  @param P: type used to gather/scatter data in and out of the message buffer
        *  @param codim: communicate entites of given codim
        *  @param if: one of the predifined interface types, throws error if it is not implemented
        *  @param level: communicate for entities on the given level
        *
        *  Implements a generic communication function sending an object of type P for each entity
        *  in the intersection of two processors. P has two methods gather and scatter that implement
        *  the protocol. Therefore P is called the "protocol class".
        */
        template<class T, template<class> class P, int codim>
        void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);
        
        /*! The new communication interface
        
        communicate objects for all codims on a given level
        */
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
        {}
        
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {}
#endif
        
        
        /** dummy collective communication */
        const CollectiveCommunication& comm () const
        {
            return ccobj_;
        }
        
        
        // **********************************************************
        // End of Interface Methods
        // **********************************************************
        
    private:

        //! compute the grid indices and ids
    void setIndices()
    {
        // //////////////////////////////////////////
        //   Create the index sets
        // //////////////////////////////////////////
        for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
            FoamGridLevelIndexSet<const FoamGrid >* p
                = new FoamGridLevelIndexSet<const FoamGrid >();
            levelIndexSets_.push_back(p);
        }
        
        for (int i=0; i<=maxLevel(); i++)
            if (levelIndexSets_[i])
                levelIndexSets_[i]->update(*this, i);
        
        //leafIndexSet_.update(*this);
#warning leafIndexSet_.update NOT called            
        
        // IdSets don't need updating

        }
         
        //! \todo Please doc me !
        CollectiveCommunication ccobj_;

    // Stores the lists of vertices, edges, elements for each level
    std::vector<tuple<std::list<FoamGridVertex>,
                      std::list<FoamGridEdge>,
                      std::list<FoamGridElement> > > entityImps_;

        //! Our set of level indices
        std::vector<FoamGridLevelIndexSet<const FoamGrid>*> levelIndexSets_;
        
        //! \todo Please doc me !
        FoamGridLeafIndexSet<const FoamGrid > leafIndexSet_;
    
        //! \todo Please doc me !
        FoamGridGlobalIdSet<const FoamGrid > globalIdSet_;
    
        //! \todo Please doc me !
        FoamGridLocalIdSet<const FoamGrid > localIdSet_;
    
}; // end Class FoamGrid




namespace Capabilities
{
    //! \todo Please doc me !
    template<int codim>
    struct hasEntity< FoamGrid, codim>
    {
        static const bool v = true;
    };
    
    
    //! \todo Please doc me !
    template <>
    struct isParallel< FoamGrid >
    {
        static const bool v = false;
    };
    
    
    //! \todo Please doc me !
    template<>
    struct hasHangingNodes< FoamGrid >
    {
        static const bool v = false;
    };

    //! \todo Please doc me !
    template<>
    struct isLevelwiseConforming< FoamGrid >
    {
        static const bool v = true;
    };

    //! \todo Please doc me !
    template<>
    struct isLeafwiseConforming< FoamGrid >
    {
        static const bool v = true;
    };
}

} // namespace Dune


#include <dune-foamgrid/foamgrid/foamgridfactory.hh>

#endif
