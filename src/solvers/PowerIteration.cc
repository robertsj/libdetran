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
                                  SP_externalsource q_e,
                                  SP_fissionsource  q_f)
  : d_input(input)
  , d_state(state)
  , d_mesh(mesh)
  , d_material(material)
  , d_quadrature(quadrature)
  , d_boundary(boundary)
  , d_fission_source(q_f)
  , d_max_iters(100)
  , d_tolerance(1e-5)
  , d_print_out(2)
  , d_print_interval(10)
  , d_aitken(false)
{
  Require(d_input);
  Require(d_state);
  Require(d_mesh);
  Require(d_material);
  Require(d_quadrature);
  Require(d_boundary);
  Require(d_fission_source);
  // Get relevant input parameters.
  if (input->check("eigen_max_iters"))
    d_max_iters = input->get<int>("eigen_max_iters");

  if (input->check("eigen_tolerance"))
    d_tolerance = input->get<double>("eigen_tolerance");

  if (input->check("eigen_print_out"))
    d_print_out = input->get<int>("eigen_print_out");

  if (input->check("eigen_print_interval"))
    d_print_interval = input->get<int>("eigen_print_interval");

  if (input->check("eigen_aitken"))
    d_aitken = input->get<int>("eigen_aitken");


  // Create multigroup solver.
  d_mg_solver = new GaussSeidel<D>(input, state, mesh, material,
                                   quadrature, boundary, q_e, q_f);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of PowerIteration.cc
//---------------------------------------------------------------------------//



