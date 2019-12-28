// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_FACTORY_HH
#define DUNE_FOAMGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for FoamGrid
    \author Oliver Sander
 */

#include <vector>
#include <memory>
#include <unordered_map>

#include <dune/common/version.hh>
#include <dune/common/function.hh>
#include <dune/common/fvector.hh>
#include <dune/common/to_unique_ptr.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

namespace Dune {

/** \brief Specialization of the generic GridFactory for FoamGrid<dimgrid, dimworld>
    */
template <int dimgrid, int dimworld, class ct>
    class GridFactoryBase
        : public GridFactoryInterface<FoamGrid<dimgrid, dimworld, ct> >
    {
    /** \brief Type used by the grid for coordinates */
    using ctype = ct;

    public:

        /** \brief Default constructor */
        GridFactoryBase()
        : factoryOwnsGrid_(true)
        {
            grid_ = new FoamGrid<dimgrid, dimworld, ctype>;
            grid_->entityImps_.resize(1);
        }

        /** \brief Constructor for a given grid object

        If you already have your grid object constructed you can
        hand it over using this constructor.

        If you construct your factory class using this constructor
        the pointer handed over to you by the method createGrid() is
        the one you supplied here.
         */
        GridFactoryBase(FoamGrid<dimgrid, dimworld, ctype>* grid)
        : factoryOwnsGrid_(false)
        {
            grid_ = grid;
            grid_->entityImps_.resize(1);
        }

        /** \brief Destructor */
        ~GridFactoryBase() override {
            if (grid_ && factoryOwnsGrid_)
                delete grid_;
        }

        /** \brief Insert a vertex into the coarse grid */
        void insertVertex(const FieldVector<ctype,dimworld>& pos) override {
            std::get<0>(grid_->entityImps_[0]).emplace_back(
              0,    // level
              pos,  // position
              grid_->getNextFreeId()
            );
            vertexArray_.push_back(&*std::get<0>(grid_->entityImps_[0]).rbegin());
        }

        /** \brief Obtain an element's insertion index
         */
        unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld, ctype>::Traits::template Codim<0>::Entity &entity) const override
        {
            return entity.impl().target_->leafIndex_;
        }

        /** \brief Obtain a vertex' insertion index
         */
        unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld, ctype>::Traits::template Codim<dimgrid>::Entity &vertex) const override
        {
            return vertex.impl().target_->leafIndex_;
        }

        /** \brief Obtain a boundary's insertion index
         */
        unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld, ctype>::LeafIntersection& intersection) const override
        {
            return intersection.boundarySegmentIndex();
        }

    protected:

        // Pointer to the grid being built
        FoamGrid<dimgrid, dimworld, ctype>* grid_;

        // True if the factory allocated the grid itself, false if the
        // grid was handed over from the outside
        bool factoryOwnsGrid_;

        /** \brief Array containing all vertices */
        std::vector<FoamGridEntityImp<0, dimgrid, dimworld, ctype>*> vertexArray_;

        /** \brief Counter that creates the boundary segment indices */
        unsigned int boundarySegmentCounter_ = 0;
    };

template <int dimgrid, int dimworld, class ct>
    class GridFactory<FoamGrid<dimgrid, dimworld, ct> >
        : public GridFactoryBase<dimgrid, dimworld, ct>
    {};

/** \brief Specialization of GridFactoryBase for 1D-FoamGrid<1, dimworld>
    */
template <int dimworld, class ct>
    class GridFactory<FoamGrid<1, dimworld, ct> >
        : public GridFactoryBase<1, dimworld, ct>
    {
        /** \brief Grid dimension */
        enum {dimgrid = 1};
        typedef ct ctype;

    public:

        GridFactory() {}

        GridFactory(FoamGrid<1, dimworld, ctype>* grid):
            GridFactoryBase<1,dimworld,ctype>(grid)
        {}

        /** \brief Insert a boundary segment.
            This is only needed if you want to control the numbering of the boundary segments
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices) override
        {
            boundarySegmentIndices_[ vertices[0] ] = this->boundarySegmentCounter_++;
        }

        /** \brief Insert a boundary segment and the boundary segment geometry
            This influences the ordering of the boundary segments.
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<BoundarySegment<dimgrid, dimworld> >& boundarySegment) override
        {
            DUNE_THROW(Dune::NotImplemented, "Parameterized boundary segments are not implemented");
        }

        /** \brief Return true if leaf intersection was inserted as boundary segment
        */
        bool wasInserted( const typename FoamGrid<dimgrid, dimworld, ctype>::LeafIntersection &intersection ) const override
        {
          if ( !intersection.boundary() || intersection.inside().level() != 0 )
            return false;

          const auto& vertex = intersection.inside().template subEntity<1>(intersection.indexInInside());
          const auto& it = boundarySegmentIndices_.find( this->insertionIndex(vertex) );
          return (it != boundarySegmentIndices_.end());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices) override {
            assert(type.isLine());
            FoamGridEntityImp<1, dimgrid, dimworld, ctype> newElement(this->vertexArray_[vertices[0]],
                                                                      this->vertexArray_[vertices[1]],
                                                                      0,
                                                                      this->grid_->getNextFreeId());

            std::get<1>(this->grid_->entityImps_[0]).push_back(newElement);

        }

        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices,
                           const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization) override
        {
            assert(type.isLine());
            FoamGridEntityImp<1, dimgrid, dimworld, ctype> newElement(this->vertexArray_[vertices[0]],
                                                                      this->vertexArray_[vertices[1]],
                                                                      0,
                                                                      this->grid_->getNextFreeId());
            // save the pointer to the element parametrization
            newElement.elementParametrization_ = elementParametrization;

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Finalize grid creation and hand over the grid

        The receiver takes responsibility of the memory allocated for the grid
        */
#if DUNE_VERSION_LT(DUNE_COMMON, 2, 7)
        FoamGrid<1, dimworld, ctype>* createGrid() override {
#else
        ToUniquePtr<FoamGrid<1, dimworld, ctype>> createGrid() override {
#endif
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (this->grid_==nullptr)
                return nullptr;

            typename std::list<FoamGridEntityImp<1, dimgrid, dimworld, ctype> >::iterator eIt    = std::get<1>(this->grid_->entityImps_[0]).begin();
            typename std::list<FoamGridEntityImp<1, dimgrid, dimworld, ctype> >::iterator eEndIt = std::get<1>(this->grid_->entityImps_[0]).end();

            for(;eIt!=eEndIt;eIt++) {

                // Get two vertices of the edge
                const FoamGridEntityImp<0, dimgrid, dimworld, ctype>* v0 = eIt->vertex_[0];
                const FoamGridEntityImp<0, dimgrid, dimworld, ctype>* v1 = eIt->vertex_[1];

                // make vertices know about edge
                // using const_cast because of the implementation of FoamGridEntityImp<1,dimgrid,dimworld>
                // the member variable vertex_ is an array with pointers to const vertices
                const_cast <FoamGridEntityImp<0, dimgrid, dimworld, ctype>*> (v0)->elements_.push_back(&*eIt);
                const_cast <FoamGridEntityImp<0, dimgrid, dimworld, ctype>*> (v1)->elements_.push_back(&*eIt);

            }

            // Create the index sets
            this->grid_->setIndices();

            // ////////////////////////////////////////////////
            //   Set the boundary ids
            // ////////////////////////////////////////////////

            // Iterate over all facets (=vertices in 1d)
            auto fIt = std::get<0>(this->grid_->entityImps_[0]).begin();
            const auto fEndIt = std::get<0>(this->grid_->entityImps_[0]).end();
            for (; fIt != fEndIt; ++fIt)
                if(fIt->elements_.size()==1) // if boundary facet
                {
                  const auto& it = boundarySegmentIndices_.find( fIt->vertex_[0]->leafIndex_ );
                  if (it != boundarySegmentIndices_.end())
                      fIt->boundarySegmentIndex_ = it->second;
                  else
                      fIt->boundarySegmentIndex_ = this->boundarySegmentCounter_++;
                }

            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimgrid, dimworld, ctype>* tmp = this->grid_;
            tmp->numBoundarySegments_ = this->boundarySegmentCounter_;
            this->grid_ = nullptr;
            return tmp;
        }

      private:
        std::unordered_map<unsigned int, unsigned int> boundarySegmentIndices_;
    };

    /** \brief Specialization of GridFactoryBase for 2D-FoamGrid<2, dimworld>
    */
template <int dimworld, class ct>
    class GridFactory<FoamGrid<2, dimworld, ct> >
        : public GridFactoryBase<2, dimworld, ct>
    {
        /** \brief Grid dimension */
        enum {dimgrid = 2};
        typedef ct ctype;

    public:

        GridFactory() {}

        GridFactory(FoamGrid<2, dimworld, ctype>* grid):
            GridFactoryBase<2, dimworld, ctype>(grid)
        {}

        /** \brief Insert a boundary segment.
            This is only needed if you want to control the numbering of the boundary segments
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices) override
        {
            std::array<unsigned int, 2> vertexIndices {{ vertices[0], vertices[1] }};

            // sort the indices
            if ( vertexIndices[0] > vertexIndices[1] )
              std::swap( vertexIndices[0], vertexIndices[1] );

            boundarySegmentIndices_[ vertexIndices ] = this->boundarySegmentCounter_++;
        }

        /** \brief Insert a boundary segment (== a line) and the boundary segment geometry
            This influences the ordering of the boundary segments.
        */
        void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<BoundarySegment<dimgrid, dimworld> >& boundarySegment) override
        {
            DUNE_THROW(Dune::NotImplemented, "Parameterized boundary segments are not implemented");
        }

        /** \brief Return true if leaf intersection was inserted as boundary segment
        */
        bool wasInserted( const typename FoamGrid<dimgrid, dimworld, ctype>::LeafIntersection &intersection ) const override
        {
          if ( !intersection.boundary() || intersection.inside().level() != 0 )
            return false;

          // Get the vertices of the intersection
          const auto refElement = ReferenceElements<ctype, dimgrid>::general(intersection.inside().type());

          const int subIdx0 = refElement.subEntity(intersection.indexInInside(), 1, /*idx*/0, dimgrid);
          const auto vertex0 = intersection.inside().template subEntity<2>( subIdx0 );
          const int subIdx1 = refElement.subEntity(intersection.indexInInside(), 1, /*idx*/1, dimgrid);
          const auto vertex1 = intersection.inside().template subEntity<2>( subIdx1 );

          std::array<unsigned int, 2> vertexIndices {{
            this->insertionIndex( vertex0 ),
            this->insertionIndex( vertex1 )
          }};

          // sort the indices
          if ( vertexIndices[0] > vertexIndices[1] )
            std::swap( vertexIndices[0], vertexIndices[1] );

          const auto& it = boundarySegmentIndices_.find( vertexIndices );
          return (it != boundarySegmentIndices_.end());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices) override {

            assert(type.isTriangle());

            FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype> newElement(0,   // level
                                       this->grid_->getNextFreeId());  // id
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        void insertElement(const GeometryType& type,
                           const std::vector<unsigned int>& vertices,
                           const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization) override
        {
            assert(type.isTriangle());
            FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype> newElement(0,   // level
                                       this->grid_->getNextFreeId());  // id
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];
            // save the pointer to the element parametrization
            newElement.elementParametrization_ = elementParametrization;

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Finalize grid creation and hand over the grid
        The receiver takes responsibility of the memory allocated for the grid
        */
#if DUNE_VERSION_LT(DUNE_COMMON, 2, 7)
        FoamGrid<dimgrid, dimworld, ctype>* createGrid() override {
#else
        ToUniquePtr<FoamGrid<dimgrid, dimworld, ctype>> createGrid() override {
#endif
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (this->grid_==nullptr)
                return nullptr;

            // ////////////////////////////////////////////////////
            //   Create the edges
            // ////////////////////////////////////////////////////

            // for convenience
            typedef FoamGridEntityImp<0, dimgrid, dimworld, ctype> FoamGridVertex;

            // For fast retrieval: a map from pairs of vertices to the edge that connects them
            std::map<std::pair<const FoamGridEntityImp<0, dimgrid, dimworld, ctype>*, const FoamGridEntityImp<0, dimgrid, dimworld, ctype>*>, FoamGridEntityImp<1, dimgrid, dimworld, ctype>*> edgeMap;

            typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype> >::iterator eIt    = std::get<dimgrid>(this->grid_->entityImps_[0]).begin();
            typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype> >::iterator eEndIt = std::get<dimgrid>(this->grid_->entityImps_[0]).end();

            for (; eIt!=eEndIt; ++eIt) {

                FoamGridEntityImp<dimgrid, dimgrid, dimworld, ctype>* element = &(*eIt);

                const auto refElement = ReferenceElements<ctype, dimgrid>::general(eIt->type());

                // Loop over all edges of this element
                for (size_t i=0; i<element->facet_.size(); ++i) {

                    // Get two vertices of the potential edge
                    const FoamGridVertex* v0 = element->vertex_[refElement.subEntity(i, 1, 0, 2)];
                    const FoamGridVertex* v1 = element->vertex_[refElement.subEntity(i, 1, 1, 2)];

                    FoamGridEntityImp<1, dimgrid, dimworld, ctype>* existingEdge = nullptr;
                    typename std::map<std::pair<const FoamGridEntityImp<0, dimgrid, dimworld, ctype>*, const FoamGridEntityImp<0, dimgrid, dimworld, ctype>*>, FoamGridEntityImp<1, dimgrid, dimworld, ctype>*>::const_iterator e = edgeMap.find(std::make_pair(v0,v1));

                    if (e != edgeMap.end()) {
                        existingEdge = e->second;
                    } else {
                        e = edgeMap.find(std::make_pair(v1,v0));
                        if (e != edgeMap.end())
                            existingEdge = e->second;
                    }

                    if (existingEdge == nullptr) {

                        // The current edge has not been inserted already.  We do that now.
                        std::get<1>(this->grid_->entityImps_[0]).push_back(FoamGridEntityImp<1, dimgrid, dimworld, ctype>(v0,
                                                                                                    v1,
                                                                                                    0, // level
                                                                                                    this->grid_->getNextFreeId() // id
                                                                                                    ));

                        existingEdge = &*std::get<1>(this->grid_->entityImps_[0]).rbegin();

                        edgeMap.insert(std::make_pair(std::make_pair(v0,v1), existingEdge));

                    }

                    // make element know about the edge
                    element->facet_[i] = existingEdge;

                    // make edge know about the element
                    existingEdge->elements_.push_back(element);
                }
            }

            // Create the index sets
            this->grid_->setIndices();


            // ////////////////////////////////////////////////
            //   Set the boundary ids
            // ////////////////////////////////////////////////

            // Iterate over all facets (=edges in 2D)
            auto fIt = std::get<1>(this->grid_->entityImps_[0]).begin();
            const auto fEndIt = std::get<1>(this->grid_->entityImps_[0]).end();
            for (; fIt!=fEndIt; ++fIt)
                if(fIt->elements_.size()==1) //if boundary facet
                {
                  std::array<unsigned int, 2> vertexIndices {{ fIt->vertex_[0]->leafIndex_, fIt->vertex_[1]->leafIndex_ }};

                  // sort the indices
                  if ( vertexIndices[0] > vertexIndices[1] )
                    std::swap( vertexIndices[0], vertexIndices[1] );

                  auto it = boundarySegmentIndices_.find( vertexIndices );
                  if (it != boundarySegmentIndices_.end()) {
                      fIt->boundarySegmentIndex_ = it->second;
                  } else { // edge was not inserted as boundary segment
                      fIt->boundarySegmentIndex_ = this->boundarySegmentCounter_++;
                  }
                }

            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimgrid, dimworld, ctype>* tmp = this->grid_;
            tmp->numBoundarySegments_ = this->boundarySegmentCounter_;
            this->grid_ = nullptr;
            return tmp;
        }

      private:
        struct HashUIntArray {
          std::size_t operator() (const std::array<unsigned int, 2>& a) const {
            return std::hash<unsigned int>{}(a[0]) ^ (std::hash<unsigned int>{}(a[1]) << 1);
          }
        };

        std::unordered_map<std::array<unsigned int, 2>, unsigned int, HashUIntArray> boundarySegmentIndices_;
    };

} // end namespace Dune

#endif
