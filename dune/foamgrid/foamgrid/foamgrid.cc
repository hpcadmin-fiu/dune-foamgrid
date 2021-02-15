#include <config.h>

#include <dune/foamgrid/foamgrid.hh>

namespace Dune {

template class FoamGrid<1,1,double>;
template class FoamGrid<1,2,double>;
template class FoamGrid<1,3,double>;
template class FoamGrid<2,2,double>;
template class FoamGrid<2,3,double>;

} // end namespace Dune
