//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  WGSolverSI.cc
 *  @brief WGSolverSI member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "WGSolverSI.hh"

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
WGSolverSI<D>::WGSolverSI(SP_state                  state,
                          SP_material               material,
                          SP_quadrature             quadrature,
                          SP_boundary               boundary,
                          const vec_externalsource &q_e,
                          SP_fissionsource          q_f,
                          bool                      multiply)
  : Base(state, material, quadrature, boundary, q_e, q_f, multiply)
{
  d_sweeper->set_update_boundary(true);
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class WGSolverSI<_1D>;
template class WGSolverSI<_2D>;
template class WGSolverSI<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of WGSolverSI.cc
//----------------------------------------------------------------------------//
