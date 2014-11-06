// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_ENTITY_HH
#define DUNE_FOAMGRID_ENTITY_HH

/** \file
* \brief The FoamGridEntity class
*/

#include <memory>

#include <dune/common/version.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/grid.hh>

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridgeometry.hh>

namespace Dune {


// Forward declarations

template<int codim, int dimentity, class GridImp>
class FoamGridEntity;

template<int codim, class GridImp>
class FoamGridEntityPointer;

template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLevelIterator;

template<class GridImp>
class FoamGridLevelIntersectionIterator;

template<class GridImp>
class FoamGridLeafIntersectionIterator;

template<class GridImp>
class FoamGridHierarchicIterator;


/** \brief The implementation of entities in a FoamGrid
*   \ingroup FoamGrid
*
*
*/
template<int codim, int dimgrid, class GridImp>
class FoamGridEntity :
    public EntityDefaultImplementation <codim, dimgrid, GridImp, FoamGridEntity>
{
    template <class GridImp_>
    friend class FoamGridLevelIndexSet;

    template <class GridImp_>
    friend class FoamGridLeafIndexSet;

    template <class GridImp_>
    friend class FoamGridLocalIdSet;

    template <class GridImp_>
    friend class FoamGridGlobalIdSet;

    friend class FoamGridEntityPointer<codim,GridImp>;


    private:

        typedef typename GridImp::ctype ctype;

        enum{dimworld = GridImp::dimensionworld};

    public:

        typedef typename GridImp::template Codim<codim>::Geometry Geometry;

        //! The type of the EntitySeed interface class
        typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;



        //! Constructor for an entity in a given grid level
    FoamGridEntity(const FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>* target) :
            target_(target)
        {}

        /** \brief Copy constructor */
        FoamGridEntity(const FoamGridEntity& original) :
            target_(original.target_)
        {}


        //! \todo Please doc me !
        FoamGridEntity& operator=(const FoamGridEntity& original)
        {
            if (this != &original)
            {
                target_ = original.target_;
            }
            return *this;
        }


        //! level of this element
        int level () const {
            return target_->level();
        }


        /** \brief The partition type for parallel computing
        */
        PartitionType partitionType () const {
            return target_->partitionType();
        }


        /** Intra-element access to entities of codimension cc > codim. Return number of entities
        * with codimension cc.
        */
        template<int cc> int count () const{
            return target_->template count<cc>();
        }


        //! geometry of this entity
        Geometry geometry () const
        {
            std::vector<FieldVector<double, dimworld> > coordinates(target_->corners());
            for (size_t i=0; i<target_->corners(); i++)
                coordinates[i] = target_->corner(i);

            return Geometry(FoamGridGeometry<dimgrid-codim, dimworld, GridImp>(target_->type(), coordinates));
        }

        //! Create EntitySeed
        EntitySeed seed () const
        {
            return EntitySeed(target_);
        }

        const FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>* target_;


        //! \todo Please doc me !
        void setToTarget(const FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>* target)
        {
            target_ = target;
        }

};




/** \brief Specialization for codim-0-entities, i.e., elements.
* \ingroup FoamGrid
*
* This class embodies the topological parts of elements of the grid.
* It has an extended interface compared to the general entity class.
* For example, Entities of codimension 0  allow to visit all neighbors.
*/
template<int dimgrid, class GridImp>
class FoamGridEntity<0, dimgrid, GridImp> :
    public EntityDefaultImplementation<0, dimgrid, GridImp, FoamGridEntity>
{

    enum {dimworld = GridImp::dimensionworld};

    public:

        typedef typename GridImp::template Codim<0>::Geometry Geometry;

        typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

        //! The Iterator over intersections on this level
        typedef FoamGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;

        //! The Iterator over intersections on the leaf level
        typedef FoamGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

        //! Iterator over descendants of the entity
        typedef FoamGridHierarchicIterator<GridImp> HierarchicIterator;

        //! The type of the EntitySeed interface class
        typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;


        //! Constructor for an entity in a given grid level
        FoamGridEntity(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* hostEntity) :
            target_(hostEntity)
        {}


        /** \brief Copy constructor */
        FoamGridEntity(const FoamGridEntity& original) :
            target_(original.target_)
        {}


        //! \todo Please doc me !
        FoamGridEntity& operator=(const FoamGridEntity& original)
        {
            if (this != &original)
            {
                target_ = original.target_;
            }
            return *this;
        }


        //! Level of this element
        int level () const
        {
            return target_->level_;
        }


        /** \brief The partition type for parallel computing */
        PartitionType partitionType () const {
            return InteriorEntity;
        }


        //! Geometry of this entity
        Geometry geometry () const
        {
            std::vector<FieldVector<double, dimworld> > coordinates(target_->corners());
            for (size_t i=0; i<target_->corners(); i++)
                coordinates[i] = target_->vertex_[i]->pos_;

            return Geometry(FoamGridGeometry<dimgrid, dimworld, GridImp>(target_->type(), coordinates));
        }

        //! Create EntitySeed
        EntitySeed seed () const
        {
            return EntitySeed(target_);
        }


        /** \brief Return the number of subEntities of codimension cc.
        */
        template<int cc>
        int count () const
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,4)
            static_assert(0<=cc && cc<=2, "Only codimensions with 0 <= cc <= 2 are valid!");
#else
            dune_static_assert(0<=cc && cc<=2, "Only codimensions with 0 <= cc <= 2 are valid!");
#endif
            if(dimgrid == 2)
                return (cc==0) ? 1 : 3;
            else 
                return (cc==0) ? 1 : 2;
        }

        /** \brief Return the number of subEntities of codimension cc.
        */
        unsigned int count (unsigned int codim) const
        {
            assert(0<=codim && codim<=2);
            if(dimgrid == 2)
                return (codim==0) ? 1 : 3;
            else 
                return (codim==0) ? 1 : 2;
        }

        /*! Return number of subentities with codimension cc.
     	 *
         * That number is (dim over (dim-codim)) times 2^codim
         */
    	unsigned int subEntities (unsigned int codim) const
    	{
      		assert(0<=codim && codim<=2);
            if(dimgrid == 2)
                return (codim==0) ? 1 : 3;
            else 
                return (codim==0) ? 1 : 2;
    	}

    /** \brief Return index of sub entity with codim = cc and local number i
     */
    int subId (int i, unsigned int codim) const {
        assert(0<=codim && codim<=dimgrid);
        if(dimgrid == 2) {
            switch (codim) {
            case 0:
                return target_->id_;
            case 1:
                return target_->edges_[i]->id_;
            case 2:
                return target_->vertex_[i]->id_;
           }
        } else {
            switch (codim) {
            case 0:
                return target_->id_;
            case 1:
                return target_->vertex_[i]->id_;
        }
        }
        DUNE_THROW(GridError, "Non-existing codimension requested!");
    }



        /** \brief Provide access to sub entity i of given codimension. Entities
        *  are numbered 0 ... count<cc>()-1
        */
        template<int codim>
        typename GridImp::template Codim<codim>::EntityPointer subEntity (int i) const{
            if (codim==0) {
                // The cast is correct when this if clause is executed
                return FoamGridEntityPointer<codim, GridImp>( (FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>*)this->target_);
            } else if (codim==1) {
                // The cast is correct when this if clause is executed
                return FoamGridEntityPointer<codim, GridImp>( (FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>*)this->target_->edges_[i]);
            } else if (codim==2) {
                // The cast is correct when this if clause is executed
                return FoamGridEntityPointer<codim, GridImp>( (FoamGridEntityImp<dimgrid-codim, dimgrid, dimworld>*)this->target_->vertex_[i]);
            }
        }


        //! First level intersection
        FoamGridLevelIntersectionIterator<GridImp> ilevelbegin () const{
            return FoamGridLevelIntersectionIterator<GridImp>(target_, 0);
        }


        //! Reference to one past the last neighbor
        FoamGridLevelIntersectionIterator<GridImp> ilevelend () const{
            return FoamGridLevelIntersectionIterator<GridImp>(target_);
        }


        //! First leaf intersection
        FoamGridLeafIntersectionIterator<GridImp> ileafbegin () const{
            if(isLeaf())
                return FoamGridLeafIntersectionIterator<GridImp>(target_,0);
            else
                return FoamGridLeafIntersectionIterator<GridImp>(target_);
        }


        //! Reference to one past the last leaf intersection
        FoamGridLeafIntersectionIterator<GridImp> ileafend () const{
            return FoamGridLeafIntersectionIterator<GridImp>(target_);
        }


        //! returns true if Entity has NO children
        bool isLeaf() const {
            return target_->isLeaf();
        }

    /** \brief Return true if this element has a father element */
    bool hasFather() const {
        return level()>0;
    }

    bool isNew() const
    {
        return target_->isNew();
    }

    bool mightVanish() const
    {
        return target_->mightVanish();
    }
        //! Inter-level access to father element on coarser grid.
        //! Assumes that meshes are nested.
        FoamGridEntityPointer<0, GridImp> father () const {
            return FoamGridEntityPointer<0, GridImp>(target_->father_);
        }


        /** \brief Location of this element relative to the reference element element of the father.
        * This is sufficient to interpolate all dofs in conforming case.
        * Nonconforming may require access to neighbors of father and
        * computations with local coordinates.
        * On the fly case is somewhat inefficient since dofs  are visited several times.
        * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
        * implementation of numerical algorithms is only done for simple discretizations.
        * Assumes that meshes are nested.
        */
        LocalGeometry geometryInFather () const {
            FoamGridEntityImp<dimgrid, dimgrid, dimworld>* father = target_->father_;
            // Check whether there really is a father
            if(father==nullptr)
                DUNE_THROW(GridError, "There is no father Element.");

            // Sanity check
            if(target_->type().isTriangle()){
                // Lookup the coordinates within the father
                // As in the refinement routine the children
                // are number as follows:
                // First come the ones located in the corner
                // ascending with the corner index.
                // Their first corner (origin in the reference simplex)
                // is always the corner that is also a corner of the father.
                // For the element with all corners on the edge midpoints of
                // the father, the corner are numbered according to the edge indices
                // of the father.
                double mapping[4][3][2] ={
                    { {0.0,0.0}, {0.5,0.0}, {0.0,0.5} },
                    { {1.0,0.0}, {0.5,0.5}, {0.5,0.0} },
                    { {0.0,1.0}, {0.0,0.5}, {0.5,0.5} },
                    { {0.5,0.0}, {0.5,0.5}, {0.0,0.5} }
                };

                std::vector<FieldVector<typename GridImp::ctype, dimgrid> >
                    coordinates(3);

                for(int corner=0; corner <3; ++corner)
                    for(int entry=0; entry <2; ++entry)
                        coordinates[corner][entry]=
                            mapping[target_->refinementIndex_][corner][entry];

                // return LocalGeomety by value
                return LocalGeometry(FoamGridGeometry<dimgrid, dimgrid, GridImp>(target_->type(),
                                                                        coordinates));
            }else{
                DUNE_THROW(NotImplemented, "geometryInFather only supported for triangles!");
            }

        }


        /** \brief Inter-level access to son elements on higher levels<=maxlevel.
        * This is provided for sparsely stored nested unstructured meshes.
        * Returns iterator to first son.
        */
        FoamGridHierarchicIterator<GridImp> hbegin (int maxLevel) const
        {
            FoamGridHierarchicIterator<GridImp> it(maxLevel);

            // Load sons of old target onto the iterator stack
            if (level()<=maxLevel && !isLeaf())
                for (size_t i=0; i<target_->nSons(); i++)
                    it.elemStack.push(target_->sons_[i]);

            GridImp::getRealImplementation(it.virtualEntity_).setToTarget((it.elemStack.empty())
                                          ? nullptr : it.elemStack.top());

            return it;
        }


        //! Returns iterator to one past the last son
        FoamGridHierarchicIterator<GridImp> hend (int maxLevel) const
        {
            return FoamGridHierarchicIterator<const GridImp>(maxLevel);
        }


        // /////////////////////////////////////////
        //   Internal stuff
        // /////////////////////////////////////////


        /** \brief Make this class point to a new FoamGridEntityImp object */
        void setToTarget(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* target)
        {
            target_ = target;
        }

        const FoamGridEntityImp<dimgrid, dimgrid ,dimworld>* target_;

    private:

        typedef typename GridImp::ctype ctype;

}; // end of FoamGridEntity codim = 0


} // namespace Dune


#endif
