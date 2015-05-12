// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_HH
#define DUNE_FOAMGRID_HH

/** \file
* \brief The FoamGrid class
*/

#include <list>
#include <set>

#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/tuples.hh>
#include <dune/common/stdstreams.hh>
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
#include "foamgrid/foamgridentityseed.hh"
#include "foamgrid/foamgridintersectioniterators.hh"
#include "foamgrid/foamgridleveliterator.hh"
#include "foamgrid/foamgridleafiterator.hh"
#include "foamgrid/foamgridhierarchiciterator.hh"
#include "foamgrid/foamgridindexsets.hh"
#include "foamgrid/foamgridviews.hh"

namespace Dune {

// Forward declaration
template <int dimgrid, int dimworld>
class FoamGrid;


/** \brief Encapsulates loads of types exported by FoamGrid */
template<int dimgrid, int dimworld>
struct FoamGridFamily
{
    typedef GridTraits<
        dimgrid,   // dimgrid
        dimworld,   // dimworld
        Dune::FoamGrid<dimgrid, dimworld>,
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
        FoamGridLevelIndexSet< const FoamGrid<dimgrid, dimworld> >,
        FoamGridLeafIndexSet< const FoamGrid<dimgrid, dimworld> >,
        FoamGridIdSet< const FoamGrid<dimgrid, dimworld> >,  // global IdSet
        unsigned int,   // global id type
        FoamGridIdSet< const FoamGrid<dimgrid, dimworld> >,  // local IdSet
        unsigned int,   // local id type
        CollectiveCommunication<Dune::FoamGrid<dimgrid, dimworld> > ,
        FoamGridLevelGridViewTraits,
        FoamGridLeafGridViewTraits,
        FoamGridEntitySeed
            > Traits;
};



/** \brief An implementation of the Dune grid interface: a 1- or 2-dimensional simplicial grid in an n-dimensional world
 *
 * \tparam dimgrid Dimension of the grid; must be either 1 or 2
 * \tparam dimworld Dimension of the world space
 */
template <int dimgrid, int dimworld>
class FoamGrid :
        public GridDefaultImplementation  <dimgrid, dimworld, double, FoamGridFamily<dimgrid, dimworld> >
{

    friend class FoamGridLevelIndexSet<const FoamGrid >;
    friend class FoamGridLeafIndexSet<const FoamGrid >;
    friend class FoamGridIdSet<const FoamGrid >;
    friend class FoamGridHierarchicIterator<const FoamGrid >;
    friend class FoamGridLevelIntersectionIterator<const FoamGrid >;
    friend class FoamGridLeafIntersectionIterator<const FoamGrid >;
    friend class FoamGridLevelIntersection<const FoamGrid >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLeafIterator;

    template <class GridType_>
    friend class GridFactory;

    template <int dimgrid_, int dimworld_>
    friend class GridFactoryBase;

    template<int codim_, int dim_, class GridImp_>
    friend class FoamGridEntity;

    template <int codim_, class GridImp_>
    friend class FoamGridEntityPointer;

    template <class GridImp_, PartitionIteratorType PiType_>
    friend class FoamGridLeafGridView;
    template <class GridImp_, PartitionIteratorType PiType_>
    friend class FoamGridLevelGridView;

    public:

    /** \brief FoamGrid is only implemented for 1 and 2 dimension */
    static_assert(dimgrid==1 || dimgrid==2, "Use FoamGrid only for 1d and 2d grids!");

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! type of the used GridFamily for this grid
    typedef FoamGridFamily<dimgrid, dimworld>  GridFamily;

    //! Exports various types belonging to this grid class
    typedef typename FoamGridFamily<dimgrid, dimworld>::Traits Traits;

    //! The type used to store coordinates
    typedef double ctype;

    /** \brief Constructor, constructs an empty grid
     */
    FoamGrid()
        : leafIndexSet_(*this),
          leafGridView_(*this),
          globalRefined(0),
          numBoundarySegments_(0)
    {
        //static_assert(dimgrid == 2, "FoamGrid currently only works for 2D in nD");
        std::fill(freeIdCounter_.begin(), freeIdCounter_.end(), 0);
    }

        //! Destructor
        ~FoamGrid()
        {
            // Delete level index sets
            for (size_t i=0; i<levelIndexSets_.size(); i++)
                if (levelIndexSets_[i])
                    delete (levelIndexSets_[i]);
        }


        //! Return maximum level defined in this grid. Levels are numbered
        //! 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const {
            return entityImps_.size()-1;
        }


        //! Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid<dimgrid, dimworld> >(Dune::get<dimgrid-codim>(entityImps_[level]).begin());
        }


        //! one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid<dimgrid, dimworld> >(Dune::get<dimgrid-codim>(entityImps_[level]).end());
        }


        //! Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid<dimgrid, dimworld> >(Dune::get<dimgrid-codim>(entityImps_[level]).begin());
        }


        //! one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid<dimgrid, dimworld> >(Dune::get<dimgrid-codim>(entityImps_[level]).end());
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
            if ((codim==2 && dimgrid==2) || (codim==1 && dimgrid==1))
                return Dune::get<0>(entityImps_[level]).size();
            if ((codim==1 && dimgrid==2))
                return Dune::get<1>(entityImps_[level]).size();
            if (codim==0)
                return Dune::get<dimgrid>(entityImps_[level]).size();

            return 0;
        }


        //! number of leaf entities per codim in this process
        int size (int codim) const{
            return leafIndexSet().size(codim);
        }


        //! number of entities per level, codim and geometry type in this process
        int size (int level, GeometryType type) const {
            return this->levelIndexSet(level).size(type);
        }


        //! number of leaf entities per codim and geometry type in this process
        int size (GeometryType type) const
        {
            return this->leafIndexSet().size(type);
        }

        /** \brief The number of boundary edges on the coarsest level */
        size_t numBoundarySegments() const
        {
            return numBoundarySegments_;
        }

        /** \brief Access to the GlobalIdSet */
        const typename Traits::GlobalIdSet& globalIdSet() const{
            return idSet_;
        }


        /** \brief Access to the LocalIdSet */
        const typename Traits::LocalIdSet& localIdSet() const{
            return idSet_;
        }


        /** \brief Access to the LevelIndexSets */
        const typename Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return *levelIndexSets_[level];
        }


        /** \brief Access to the LeafIndexSet */
        const typename Traits::LeafIndexSet& leafIndexSet() const
        {
            return leafIndexSet_;
        }


        //! View for the leaf grid
        template<PartitionIteratorType pitype>
        typename Traits::template Partition<pitype>::LeafGridView
        leafGridView() const {
          typedef typename Traits::template Partition<pitype>::LeafGridView View;
          return View(leafGridView_);
        }

        //! View for the leaf grid for All_Partition
        typename Traits::template Partition<All_Partition>::LeafGridView
        leafGridView() const
        {
          typedef typename Traits::template Partition<All_Partition>::LeafGridView View;
          return View(leafGridView_);
        }


        /** \brief Create EntityPointer from EnitySeed */
        template < class EntitySeed >
        static typename Traits::template Codim<EntitySeed::codimension>::EntityPointer
        entityPointer(const EntitySeed& seed)
        {
            typedef FoamGridEntityPointer<EntitySeed::codimension, const FoamGrid> EntityPointerImpl;
            typedef typename Traits::template Codim<EntitySeed::codimension>::EntityPointer EntityPointer;

            return EntityPointer(EntityPointerImpl(FoamGrid::getRealImplementation(seed).getImplementationPointer()));
        }

        /** \brief Create an Entity from an EntitySeed */
        template <class EntitySeed>
        static typename Traits::template Codim<EntitySeed::codimension>::Entity
        entity(const EntitySeed& seed)
        {
          const int codim = EntitySeed::codimension;
          typedef typename Traits::template Codim<codim>::Entity Entity;
          return Entity(FoamGridEntity<codim, dimgrid, const FoamGrid>(FoamGrid::getRealImplementation(seed).getImplementationPointer()));
        }


        /** @name Grid Refinement Methods */
        /*@{*/


        /** \brief Refine the grid uniformly
         * \param refCount Number of times the grid is to be refined uniformly
        */
        void globalRefine (int refCount = 1);

        /** \brief Mark entity for refinement
        *
        * This only works for entities of codim 0.
        * The parameter is currently ignored
        *
        * \return <ul>
        * <li> true, if marking was successful </li>
        * <li> false, if marking was not possible </li>
        * </ul>
        */
        bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e)
        {
            if (!e.isLeaf())
                return false;

            /** \todo Why do I need those const_casts here? */
            if (refCount>=1)
                const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(e).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::REFINE;
            else if (refCount<0)
                const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(e).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::COARSEN;
            else
                const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(e).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;

            return true;
        }

        /** \brief Return refinement mark for entity
        *
        * \return refinement mark (1,0,-1)
        */
        int getMark(const typename Traits::template Codim<0>::Entity & e) const
        {
            if (this->getRealImplementation(e).target_->markState_ == FoamGridEntityImp<dimgrid, dimgrid, dimworld>::REFINE)
                return 1;
            if (this->getRealImplementation(e).target_->markState_ == FoamGridEntityImp<dimgrid, dimgrid, dimworld>::COARSEN)
                return -1;

            return 0;
        }

        //! \brief Book-keeping routine to be called before adaptation
        bool preAdapt();

        //! Triggers the grid refinement process
        bool adapt();

        /** \brief Clean up refinement markers */
        void postAdapt();

        /*@}*/

        /** @name Methods for parallel computations */
        /*@{*/

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
        const typename Traits::CollectiveCommunication& comm () const
        {
            return ccobj_;
        }
        /*@}*/


        // **********************************************************
        // End of Interface Methods
        // **********************************************************

    private:

    //! \brief erases pointers in father elements to vanished entities of the element
    void erasePointersToEntities(std::list<FoamGridEntityImp<dimgrid, dimgrid ,dimworld> >& elements);

    //! \brief Erase Entities from memory that vanished due to coarsening.
    //! \warning This method has to be called first for i=0.
    //! \tparam i The dimension of the entities.
    //! \param  levelEntities The vector with the level entitied
    template<int i>
    void eraseVanishedEntities(std::list<FoamGridEntityImp<i, dimgrid, dimworld> >& levelEntities);

    //! \brief Coarsen an Element
    //! \param element The element to coarsen
    void coarsenSimplexElement(FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element);

    //! \brief refine an Element
    //! \param element The element to refine
    //! \param refCount How many times to refine the element
    void refineSimplexElement(FoamGridEntityImp<2, dimgrid, dimworld>& element,
                       int refCount);
    //! Overloaded function for the 1d case
    void refineSimplexElement(FoamGridEntityImp<1, dimgrid, dimworld>& element,
                       int refCount);

    /**
     * \brief Overwrites the neighbours of this and descendant edges
     *
     * After returning all neighbours the previously pointed to the
     * father will point to the son element
     * \param edge The edge to start overwriting with.
     * \param son The son element to substitute the father with.
     * \param father Pointer to the father element that is to be substituted.
     */
    void overwriteFineLevelNeighbours(FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>& edge,
                                      FoamGridEntityImp<dimgrid, dimgrid, dimworld>* son,
                                      FoamGridEntityImp<dimgrid, dimgrid, dimworld>* father);

    //! compute the grid indices and ids
    void setIndices();

    //! Collective communication interface
    typename Traits::CollectiveCommunication ccobj_;

    // Stores the lists of vertices, edges, elements for each level
    // std::vector<tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
    // std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >,
    // std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> > > > entityImps_;

    // conditional typename depending on dimension of grid (1 or 2)
    typedef typename std::conditional<
                      dimgrid==2,
                      typename std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<2, dimgrid, dimworld> > > >,
                      typename std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> > > >
                   >::type EntityImps;

    EntityImps entityImps_;

    typedef typename std::conditional<
                        dimgrid==2,
                        std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                     std::list<FoamGridEntityImp<1, dimgrid, dimworld> >,
                                     std::list<FoamGridEntityImp<2, dimgrid, dimworld> > >,
                        std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                     std::list<FoamGridEntityImp<1, dimgrid, dimworld> > >
                   >::type EntityTuple;

    //! Our set of level indices
    std::vector<FoamGridLevelIndexSet<const FoamGrid>*> levelIndexSets_;

    //! The leaf index set
    FoamGridLeafIndexSet<const FoamGrid > leafIndexSet_;

    // The leaf grid view
    FoamGridLeafGridView<const FoamGrid, All_Partition> leafGridView_;


    //! The id set
    FoamGridIdSet<const FoamGrid > idSet_;

    /** \brief Counters that always provide the next free id for each dimension */
    array<unsigned int, dimgrid+1> freeIdCounter_;

    /** \brief How many times was the leaf level globally refined. */
    int globalRefined;

    /** \brief The number of boundary segements. */
    std::size_t numBoundarySegments_;

    // True if the last call to preadapt returned true
    bool willCoarsen;
}; // end Class FoamGrid

#include <dune/foamgrid/foamgrid/foamgrid.cc>


namespace Capabilities
{
    /** \brief True if the grid implements entities of a given codim.
      *
      * FoamGrid implements all codimensions, hence this is always true
      */
    template<int dimgrid, int dimworld, int codim>
    struct hasEntity< FoamGrid<dimgrid, dimworld>, codim>
    {
        static const bool v = true;
    };


    /** \brief True if the grid can be run on a distributed machine
      */
    template <int dimgrid, int dimworld>
    struct isParallel< FoamGrid<dimgrid, dimworld> >
    {
        static const bool v = false;
    };


    //! \todo Please doc me !
    template <int dimgrid, int dimworld>
    struct isLevelwiseConforming< FoamGrid<dimgrid, dimworld> >
    {
        static const bool v = false;
    };

    //! \todo Please doc me !
    template <int dimgrid, int dimworld>
    struct isLeafwiseConforming< FoamGrid<dimgrid, dimworld> >
    {
        static const bool v = false;
    };
}

} // namespace Dune


// The factory should be present whenever the user includes foamgrid.hh.
// However since the factory needs to know the grid the include directive
// comes here at the end.
#include <dune/foamgrid/foamgrid/foamgridfactory.hh>

#endif
