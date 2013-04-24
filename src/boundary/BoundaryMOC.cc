//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryMOC.cc
 *  @brief  BoundaryMOC member definitions.
 *  @author Jeremy Roberts
 *  @date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#include "boundary/BoundaryMOC.t.hh"
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

BOUNDARY_TEMPLATE class BOUNDARY_EXPORT BoundaryMOC<_1D>;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT BoundaryMOC<_2D>;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT BoundaryMOC<_3D>;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT detran_utilities::SP<BoundaryMOC<_1D> >;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT detran_utilities::SP<BoundaryMOC<_2D> >;
BOUNDARY_TEMPLATE class BOUNDARY_EXPORT detran_utilities::SP<BoundaryMOC<_3D> >;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.cc
//---------------------------------------------------------------------------//
