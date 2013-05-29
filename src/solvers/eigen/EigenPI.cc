//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenPI.cc
 *  @brief EigenPI class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenPI.hh"
#include <iostream>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
EigenPI<D>::EigenPI(SP_mg_solver mg_solver)
  : Base(mg_solver)
  , d_aitken(false)
  , d_omega(1.0)
{
  if (d_input->check("eigen_pi_aitken"))
    d_aitken = d_input->template get<int>("eigen_pi_aitken");

  if (d_input->check("eigen_pi_omega"))
    d_omega = d_input->template get<double>("eigen_pi_omega");
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class EigenPI<_1D>;
template class EigenPI<_2D>;
template class EigenPI<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of EigenPI.cc
//----------------------------------------------------------------------------//
