//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PowerIteration.cc
 *  @author robertsj
 *  @date   Apr 10, 2012
 *  @brief  PowerIteration class definition.
 */
//---------------------------------------------------------------------------//

// Detran
#include "PowerIteration.hh"

// System
#include <iostream>

namespace detran
{

// Constructor
template <class D>
PowerIteration<D>::PowerIteration(SP_mg_solver mg_solver)
  : Base(mg_solver)
  , d_aitken(false)
{

  if (d_input->check("eigen_aitken"))
    d_aitken = d_input->template get<int>("eigen_aitken");

}

template class PowerIteration<_1D>;
template class PowerIteration<_2D>;
template class PowerIteration<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of PowerIteration.cc
//---------------------------------------------------------------------------//



