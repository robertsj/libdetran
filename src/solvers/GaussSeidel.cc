//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussSeidel.cc
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  GaussSeidel member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "GaussSeidel.hh"

// System
#include <iostream>
#include <string>

namespace detran
{

// Constructor
template <class D>
GaussSeidel<D>::GaussSeidel(SP_input          input,
                            SP_state          state,
                            SP_mesh           mesh,
                            SP_material       material,
                            SP_quadrature     quadrature,
                            SP_boundary       boundary,
                            SP_externalsource q_e,
                            SP_fissionsource  q_f)
  : Base(input, state, mesh, material, quadrature, boundary, q_e, q_f)
  , d_norm_type("Linf")
{
  if (d_input->check("outer_norm_type"))
    d_norm_type = input->template get<std::string>("outer_norm_type");

  // Post conditions
  Ensure(d_norm_type == "Linf" or d_norm_type == "L1" or d_norm_type == "L2");
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of GaussSeidel.cc
//---------------------------------------------------------------------------//

