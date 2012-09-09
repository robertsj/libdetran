//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIteration.cc
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  PowerIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
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
PowerIteration<D>::PowerIteration(SP_input          input,
                                  SP_state          state,
                                  SP_mesh           mesh,
                                  SP_material       material,
                                  SP_quadrature     quadrature,
                                  SP_boundary       boundary,
                                  SP_fissionsource  q_f)
  : Base(input, state, mesh, material, quadrature, boundary, q_f)
  , d_aitken(false)
{

  if (input->check("eigen_aitken"))
    d_aitken = input->template get<int>("eigen_aitken");

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of PowerIteration.cc
//---------------------------------------------------------------------------//



