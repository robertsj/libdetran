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

template class BoundaryFactory<_1D, BoundaryDiffusion>;
template class BoundaryFactory<_2D, BoundaryDiffusion>;
template class BoundaryFactory<_3D, BoundaryDiffusion>;
template class BoundaryFactory<_1D, BoundarySN>;
template class BoundaryFactory<_2D, BoundarySN>;
template class BoundaryFactory<_3D, BoundarySN>;
template class BoundaryFactory<_1D, BoundaryMOC>;
template class BoundaryFactory<_2D, BoundaryMOC>;
template class BoundaryFactory<_3D, BoundaryMOC>;

} // end namespace detran


