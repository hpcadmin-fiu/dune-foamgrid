// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_FACTORY_HH
#define DUNE_FOAMGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for FoamGrid
    \author Oliver Sander
 */

#include <vector>
#include <map>

#include <dune/common/fvector.hh>

#include <dune/grid/common/gridfactory.hh>

#include "../foamgrid.hh"

namespace Dune {

    /** \brief Specialization of the generic GridFactory for FoamGrid

    */
    template <int dimworld>
    class GridFactory<FoamGrid<dimworld> >
        : public GridFactoryInterface<FoamGrid<dimworld> > {

        /** \brief Type used by the grid for coordinates */
        typedef typename FoamGrid<dimworld>::ctype ctype;

        typedef typename std::map<FieldVector<ctype,1>, unsigned int>::iterator VertexIterator;

        enum {dim = FoamGrid<dimworld>::dimension};

    public:

        /** \brief Default constructor */
        GridFactory()
            : factoryOwnsGrid_(true),
              vertexIndex_(0)
        {
            grid_ = new FoamGrid<dimworld>;

            grid_->entityImps_.resize(1);
        }

        /** \brief Constructor for a given grid object

        If you already have your grid object constructed you can
        hand it over using this constructor.

        If you construct your factory class using this constructor
        the pointer handed over to you by the method createGrid() is
        the one you supplied here.
         */
        GridFactory(FoamGrid<dimworld>* grid)
            : factoryOwnsGrid_(false),
              vertexIndex_(0)
        {
            grid_ = grid;

            grid_->entityImps_.resize(1);
        }

        /** \brief Destructor */
        virtual ~GridFactory() {
            if (grid_ && factoryOwnsGrid_)
                delete grid_;
        }

        /** \brief Insert a vertex into the coarse grid */
        virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) {
            Dune::get<0>(grid_->entityImps_[0]).push_back(FoamGridEntityImp<0,dimworld> (0,   // level
                                                                         pos,  // position
                                                                         grid_->freeIdCounter_[0]++));
            vertexArray_.push_back(&*Dune::get<0>(grid_->entityImps_[0]).rbegin());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices) {
	    assert(type.isLine());
 	    FoamGridEntityImp<1,dimworld> newElement(vertexArray_[vertices[0]],vertexArray_[vertices[1]],0,grid_->freeIdCounter_[1]++);
	    Dune::get<1>(grid_->entityImps_[0]).push_back(newElement);

        }

        /** \brief Insert a boundary segment.

        This is only needed if you want to control the numbering of the boundary segments
        */
        virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices) {
            DUNE_THROW(Dune::NotImplemented, "insertBoundarySegment not implemented yet!");
        }

        /** \brief Insert a boundary segment (== a line) and the boundary segment geometry
         *
            This influences the ordering of the boundary segments.
            Currently, the BoundarySegment object does not actually have any effect.
        */
        virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                           const shared_ptr<BoundarySegment<2, dimworld> > boundarySegment)
        {
            insertBoundarySegment(vertices);
        }

        /** \brief Finalize grid creation and hand over the grid

        The receiver takes responsibility of the memory allocated for the grid
        */
        virtual FoamGrid<dimworld>* createGrid() {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (grid_==nullptr)
                return nullptr;

	        typename std::list<FoamGridEntityImp<1,dimworld> >::iterator eIt    = Dune::get<1>(grid_->entityImps_[0]).begin();
		typename std::list<FoamGridEntityImp<1,dimworld> >::iterator eEndIt = Dune::get<1>(grid_->entityImps_[0]).end();

		for(;eIt!=eEndIt;eIt++) {

                    // Get two vertices of the edge
                    const FoamGridEntityImp<0,dimworld>* v0 = eIt->vertex_[0];
                    const FoamGridEntityImp<0,dimworld>* v1 = eIt->vertex_[1];

		    // make vertices know about edge
		    // using const_cast because of the implementation of FoamGridEntityImp<1,dimworld>
		    // the member variable vertex_ is an array with pointers to const vertices
		    const_cast <FoamGridEntityImp<0,dimworld>*> (v0)->elements_.push_back(&*eIt);
		    const_cast <FoamGridEntityImp<0,dimworld>*> (v1)->elements_.push_back(&*eIt);

                }

            // Create the index sets
            grid_->setIndices();


            // ////////////////////////////////////////////////
            //   Set the boundary ids
            //  \todo It must be possible to set them by hand
            // ////////////////////////////////////////////////

            unsigned int boundaryIdCounter = 0;

            for (typename std::list<FoamGridEntityImp<0,dimworld> >::iterator it = Dune::get<0>(grid_->entityImps_[0]).begin();
                 it != Dune::get<0>(grid_->entityImps_[0]).end();
                 ++it)
                if(it->elements_.size()==1)
                    it->boundaryId_ = boundaryIdCounter++;

            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimworld>* tmp = grid_;
            tmp->numBoundarySegments_ = boundaryIdCounter;
            grid_ = nullptr;
            return tmp;
        }

    private:

        // Initialize the grid structure in UG
        void createBegin();

        // Pointer to the grid being built
        FoamGrid<dimworld>* grid_;

        // True if the factory allocated the grid itself, false if the
        // grid was handed over from the outside
        bool factoryOwnsGrid_;

        /** \brief Counter that creates the vertex indices */
        unsigned int vertexIndex_;

        std::vector<FoamGridEntityImp<0,dimworld>*> vertexArray_;

    };

}

#endif
