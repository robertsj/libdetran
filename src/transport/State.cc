//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   State.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  State member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "State.hh"

namespace detran
{

// Constructor.
State::State(SP_input        input,
             SP_mesh         mesh,
             SP_quadrature   quadrature)
  : d_input(input)
  , d_mesh(mesh)
  , d_quadrature(quadrature)
  , d_number_groups(input->get<int>("number_groups"))
  , d_store_angular_flux(false)
  , d_moments(d_number_groups, vec_dbl(mesh->number_cells(), 0.0))
  , d_eigenvalue(0.0)
{
  Require(input);
  Require(mesh);
  Require(quadrature);
  Require(d_number_groups > 0);

  // Allocate angular flux vectors if needed.
  int store_psi;
  if (input->check("store_angular_flux"))
    store_psi = input->get<int>("store_angular_flux");
  if (store_psi > 0)
  {
    d_store_angular_flux = true;
    d_angular_flux.resize(d_number_groups,
                          vec2_dbl(quadrature->number_angles(),
                                   vec_dbl(mesh->number_cells(),
                                           0.0)));
  }

}


}

//---------------------------------------------------------------------------//
//              end of State.cc
//---------------------------------------------------------------------------//
