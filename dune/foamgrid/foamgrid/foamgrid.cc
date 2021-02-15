// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/foamgrid/foamgrid.hh>

namespace Dune {

template class FoamGrid<1, 1, double>;
template class FoamGrid<1, 2, double>;
template class FoamGrid<1, 3, double>;
template class FoamGrid<2, 2, double>;
template class FoamGrid<2, 3, double>;

} // end namespace Dune
