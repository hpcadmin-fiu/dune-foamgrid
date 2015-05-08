// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FOAMGRID_NULL_ITERATORS_HH
#define DUNE_FOAMGRID_NULL_ITERATORS_HH

namespace Dune {

  // forward declaration of the entity implementation
  template <int dimentity, int dimgrid, int dimworld> class FoamGridEntityImp;

  template <int dim, int dimgrid, int dimworld>
  class FoamGridNullIteratorFactory {};

  template <int dimgrid, int dimworld>
  class FoamGridNullIteratorFactory<0, dimgrid, dimworld> 
  {
  public:
    static typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> >::iterator null() 
    {
      return emptyList_.end();
    }

  private:
    static typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> > emptyList_;
  };

  template <int dimgrid, int dimworld>
  class FoamGridNullIteratorFactory<1, dimgrid, dimworld> 
  {
  public:
    static typename std::list<FoamGridEntityImp<1, dimgrid, dimworld> >::iterator null() 
    {
      return emptyList_.end();
    }

  private:
    static typename std::list<FoamGridEntityImp<1, dimgrid, dimworld> > emptyList_;
  };

  template <int dimgrid, int dimworld>
  class FoamGridNullIteratorFactory<2, dimgrid, dimworld> 
  {
  public:
    static typename std::list<FoamGridEntityImp<2, dimgrid, dimworld> >::iterator null() 
    {
      return emptyList_.end();
    }

  private:
    static typename std::list<FoamGridEntityImp<2, dimgrid, dimworld> > emptyList_;
  };
} // end namespace Dune 

template<int dimgrid, int dimworld>
std::list<Dune::FoamGridEntityImp<0, dimgrid, dimworld> > Dune::FoamGridNullIteratorFactory<0, dimgrid, dimworld>::emptyList_;
template<int dimgrid, int dimworld>
std::list<Dune::FoamGridEntityImp<1, dimgrid, dimworld> > Dune::FoamGridNullIteratorFactory<1, dimgrid, dimworld>::emptyList_;
template<int dimgrid, int dimworld>
std::list<Dune::FoamGridEntityImp<2, dimgrid, dimworld> > Dune::FoamGridNullIteratorFactory<2, dimgrid, dimworld>::emptyList_;

#endif
