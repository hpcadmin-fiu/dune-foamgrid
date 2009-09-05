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
#include <dune-foamgrid/foamgrid.hh>

namespace Dune {

    /** \brief Specialization of the generic GridFactory for FoamGrid
        
    */
    template <>
    class GridFactory<FoamGrid> 
        : public GridFactoryInterface<FoamGrid> {

        /** \brief Type used by the grid for coordinates */
        typedef FoamGrid::ctype ctype;

        typedef std::map<FieldVector<ctype,1>, unsigned int >::iterator VertexIterator;
        
        enum {dim = FoamGrid::dimension};

    public:

        /** \brief Default constructor */
        GridFactory()
            : factoryOwnsGrid_(true),
              vertexIndex_(0)
        {
            grid_ = new FoamGrid;

            grid_->entityImps_.resize(1);
        }

        /** \brief Constructor for a given grid object 

        If you already have your grid object constructed you can
        hand it over using this constructor.

        If you construct your factory class using this constructor
        the pointer handed over to you by the method createGrid() is
        the one you supplied here.
         */
        GridFactory(FoamGrid* grid)
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
        virtual void insertVertex(const FieldVector<ctype,3>& pos) {
            Dune::get<0>(grid_->entityImps_[0]).push_back(FoamGridVertex(0,pos));
            vertexArray_.push_back(&*Dune::get<0>(grid_->entityImps_[0]).rbegin());
        }

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices) {

            assert(type.isTriangle());

            FoamGridElement newElement(0,0);
            newElement.vertex_[0] = vertexArray_[vertices[0]];
            newElement.vertex_[1] = vertexArray_[vertices[1]];
            newElement.vertex_[2] = vertexArray_[vertices[2]];

            Dune::get<dim>(grid_->entityImps_[0]).push_back(newElement);

            
        }

        /** \brief Finalize grid creation and hand over the grid

        The receiver takes responsibility of the memory allocated for the grid
        */
        virtual FoamGrid* createGrid() {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (grid_==NULL)
                return NULL;

            // Create the index sets
            grid_->setIndices();
            
            Dune::FoamGrid* tmp = grid_;
            grid_ = NULL;
            return tmp;
        }

    private:

        // Initialize the grid structure in UG
        void createBegin();

        // Pointer to the grid being built
        FoamGrid* grid_;

        // True if the factory allocated the grid itself, false if the
        // grid was handed over from the outside
        bool factoryOwnsGrid_;

        /** \brief Counter that creates the vertex indices */
        unsigned int vertexIndex_;

        std::vector<FoamGridVertex*> vertexArray_;

    };

}

#endif
