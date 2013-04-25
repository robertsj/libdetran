//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryFactory.cc
 *  @author robertsj
 *  @date   Jan 30, 2013
 *  @brief  BoundaryFactory class definition.
 */
//---------------------------------------------------------------------------//

#include "BoundaryFactory.t.hh"

namespace detran
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_1D, BoundaryDiffusion>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_2D, BoundaryDiffusion>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_3D, BoundaryDiffusion>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_1D, BoundarySN>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_2D, BoundarySN>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_3D, BoundarySN>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_1D, BoundaryMOC>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_2D, BoundaryMOC>)
BOUNDARY_TEMPLATE_EXPORT(BoundaryFactory<_3D, BoundaryMOC>)

} // end namespace detran


