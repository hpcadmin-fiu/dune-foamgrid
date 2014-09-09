// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FOAMGRID_FOAMGRID_FOAMGRIDVIEWS_HH
#define DUNE_FOAMGRID_FOAMGRID_FOAMGRIDVIEWS_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include "foamgridindexsets.hh"

namespace Dune
{

    template< class GridImp, PartitionIteratorType pitype >
    class FoamGridLevelGridView;

    template< class GridImp, PartitionIteratorType pitype >
    class FoamGridLeafGridView;


    /** \brief Collect several types associated to FoamGrid LevelGridViews */
    template< class GridImp, PartitionIteratorType pitype >
    struct FoamGridLevelGridViewTraits
    {
      typedef FoamGridLevelGridView< GridImp, pitype > GridViewImp;

      /** \brief type of the grid */
      typedef typename remove_const<GridImp>::type Grid;

      /** \brief type of the index set */
      typedef typename Grid :: Traits :: LevelIndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Grid :: Traits :: LevelIntersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Grid :: Traits :: LevelIntersectionIterator
      IntersectionIterator;

      /** \brief type of the collective communication */
      typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

      template< int cd >
      struct Codim
      {
        typedef typename Grid :: Traits
        :: template Codim< cd > :: template Partition< pitype > :: LevelIterator
        Iterator;

        typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;
        typedef typename Grid :: Traits :: template Codim< cd > :: EntityPointer
        EntityPointer;

        typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
        typedef typename Grid :: template Codim< cd > :: LocalGeometry
        LocalGeometry;

        /** \brief Define types needed to iterate over entities of a given partition type */
        template< PartitionIteratorType pit >
        struct Partition
        {
          /** \brief iterator over a given codim and partition type */
          typedef typename Grid :: template Codim< cd >
          :: template Partition< pit > :: LevelIterator
          Iterator;
        };
      };

      enum { conforming = Capabilities :: isLevelwiseConforming< Grid > :: v };
    };


    /** \brief Implementation class of LevelGridViews for FoamGrid */
    template< class GridImp, PartitionIteratorType pitype >
    class FoamGridLevelGridView
    {
    public:
      typedef FoamGridLevelGridViewTraits<GridImp,pitype> Traits;

      /** \brief type of the grid */
      typedef typename Traits::Grid Grid;

      /** \brief type of the index set */
      typedef typename Traits :: IndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Traits :: Intersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Traits :: IntersectionIterator IntersectionIterator;

      /** \brief type of the collective communication */
      typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

      /** \brief Codim Structure */
      template< int cd >
      struct Codim : public Traits :: template Codim<cd> {};

      enum { conforming = Traits :: conforming };

      FoamGridLevelGridView ( const Grid &grid, int level )
      : grid_( &grid ),
      level_( level )
      {}

      /** \brief obtain a const reference to the underlying hierarchic grid */
      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      /** \brief obtain the index set */
      const IndexSet &indexSet () const
      {
        return grid().levelIndexSet( level_ );
      }

      /** \brief obtain number of entities in a given codimension */
      int size ( int codim ) const
      {
        return grid().size( level_, codim );
      }

      /** \brief obtain number of entities with a given geometry type */
      int size ( const GeometryType &type ) const
      {
        return grid().size( level_, type );
      }

      /** \brief obtain begin iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator begin () const
      {
        return grid().template lbegin< cd, pitype >( level_ );
      }

      /** \brief obtain begin iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
      {
        return grid().template lbegin< cd, pit >( level_ );
      }

      /** \brief obtain end iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator end () const
      {
        return grid().template lend< cd, pitype >( level_ );
      }

      /** \brief obtain end iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator end () const
      {
        return grid().template lend< cd, pit >( level_ );
      }

      /** \brief obtain begin intersection iterator with respect to this view */
      IntersectionIterator
      ibegin ( const typename Codim< 0 > :: Entity &entity ) const
      {
        return grid().getRealImplementation(entity).ilevelbegin();
      }

      /** \brief obtain end intersection iterator with respect to this view */
      IntersectionIterator
      iend ( const typename Codim< 0 > :: Entity &entity ) const
      {
        //return entity.impl().ilevelend();
        return grid().getRealImplementation(entity).ilevelend();
      }

      /** \brief obtain collective communication object */
      const CollectiveCommunication &comm () const
      {
        return grid().comm();
      }

      /** \brief Return size of the overlap region for a given codim on the grid view.  */
      int overlapSize(int codim) const
      {
        return grid().overlapSize(level_, codim);
      }

      /** \brief Return size of the ghost region for a given codim on the grid view.  */
      int ghostSize(int codim) const
      {
        return grid().ghostSize(level_, codim);
      }

      /** communicate data on this view */
      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                         InterfaceType iftype,
                         CommunicationDirection dir ) const
                         {
                           return grid().communicate( data, iftype, dir, level_ );
                         }

    private:
      const Grid *grid_;
      int level_;
    };


    template< class GridImp, PartitionIteratorType pitype >
    struct FoamGridLeafGridViewTraits {
      typedef FoamGridLeafGridView< GridImp, pitype > GridViewImp;

      /** \brief type of the grid */
      typedef typename remove_const<GridImp>::type Grid;

      /** \brief type of the index set */
      typedef typename Grid :: Traits :: LeafIndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Grid :: Traits :: LeafIntersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Grid :: Traits :: LeafIntersectionIterator
      IntersectionIterator;

      /** \brief type of the collective communication */
      typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

      template< int cd >
      struct Codim
      {
        typedef typename Grid :: Traits
        :: template Codim< cd > :: template Partition< pitype > :: LeafIterator
        Iterator;

        typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;
        typedef typename Grid :: Traits :: template Codim< cd > :: EntityPointer
        EntityPointer;

        typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
        typedef typename Grid :: template Codim< cd > :: LocalGeometry
        LocalGeometry;

        /** \brief Define types needed to iterate over entities of a given partition type */
        template <PartitionIteratorType pit >
        struct Partition
        {
          /** \brief iterator over a given codim and partition type */
          typedef typename Grid :: template Codim< cd >
          :: template Partition< pit > :: LeafIterator
          Iterator;
        };
      };

      enum { conforming = Capabilities :: isLeafwiseConforming< Grid > :: v };
    };


    /** \brief Implementation class of LevelGridViews for FoamGrid */
    template< class GridImp, PartitionIteratorType pitype >
    class FoamGridLeafGridView
    {
      // Grid must be able to call the 'update' method of the index set
      friend GridImp;

    public:
      typedef FoamGridLeafGridViewTraits<GridImp,pitype> Traits;

      /** \brief type of the grid */
      typedef typename Traits::Grid Grid;

      /** \brief type of the index set */
      typedef typename Traits :: IndexSet IndexSet;

      /** \brief type of the intersection */
      typedef typename Traits :: Intersection Intersection;

      /** \brief type of the intersection iterator */
      typedef typename Traits :: IntersectionIterator IntersectionIterator;

      /** \brief type of the collective communication */
      typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

      /** \brief Codim Structure */
      template< int cd >
      struct Codim : public Traits :: template Codim<cd> {};

      enum { conforming = Traits :: conforming };

    public:
      FoamGridLeafGridView ( const Grid &grid )
      : grid_( &grid )
      {}

      /** \brief Copy constructor */
      FoamGridLeafGridView ( const FoamGridLeafGridView<GridImp,pitype>& other )
      : grid_( other.grid_ ),
        indexSet_(other.indexSet_)
      {}

      /** \brief obtain a const reference to the underlying hierarchic grid */
      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      /** \brief obtain the index set */
      const IndexSet &indexSet () const
      {
        return indexSet_;
      }

      /** \brief obtain number of entities in a given codimension */
      int size ( int codim ) const
      {
        return indexSet_.size( codim );
      }

      /** \brief obtain number of entities with a given geometry type */
      int size ( const GeometryType &type ) const
      {
        return indexSet_.size( type );
      }

      /** \brief obtain begin iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator begin () const
      {
        return grid().template leafbegin< cd, pitype >();
      }

      /** \brief obtain begin iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
      {
        return grid().template leafbegin< cd, pit >();
      }

      /** \brief obtain end iterator for this view */
      template< int cd >
      typename Codim< cd > :: Iterator end () const
      {
        return grid().template leafend< cd, pitype >();
      }

      /** \brief obtain end iterator for this view */
      template< int cd, PartitionIteratorType pit >
      typename Codim< cd > :: template Partition< pit > :: Iterator end () const
      {
        return grid().template leafend< cd, pit >();
      }

      /** \brief obtain begin intersection iterator with respect to this view */
      IntersectionIterator
      ibegin ( const typename Codim< 0 > :: Entity &entity ) const
      {
        return grid().getRealImplementation(entity).ileafbegin();
      }

      /** \brief obtain end intersection iterator with respect to this view */
      IntersectionIterator
      iend ( const typename Codim< 0 > :: Entity &entity ) const
      {
        return grid().getRealImplementation(entity).ileafend();
      }

      /** \brief obtain collective communication object */
      const CollectiveCommunication &comm () const
      {
        return grid().comm();
      }

      /** \brief Return size of the overlap region for a given codim on the grid view.  */
      int overlapSize(int codim) const
      {
        return 0;
      }

      /** \brief Return size of the ghost region for a given codim on the grid view.  */
      int ghostSize(int codim) const
      {
        return 0;
      }

      /** communicate data on this view */
      template< class DataHandleImp, class DataType >
      void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                         InterfaceType iftype,
                         CommunicationDirection dir ) const
                         {
                           return grid().communicate( data, iftype, dir );
                         }

    private:
      const Grid *grid_;
      FoamGridLeafIndexSet<GridImp> indexSet_;
    };

}

#endif // #ifndef DUNE_FOAMGRID_FOAMGRID_FOAMGRIDVIEWS_HH
