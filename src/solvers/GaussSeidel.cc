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
#include "SourceIteration.hh"

// System
#include <iostream>

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
  : d_input(input)
  , d_state(state)
  , d_mesh(mesh)
  , d_material(material)
  , d_quadrature(quadrature)
  , d_boundary(boundary)
  , d_downscatter(false)
  , d_max_iters(100)
  , d_tolerance(1e-5)
  , d_print_out(2)
  , d_print_interval(10)
{
  Require(d_input);
  Require(d_state);
  Require(d_mesh);
  Require(d_material);
  Require(d_quadrature);
  Require(d_boundary);

  Assert(d_input->check("number_groups"));
  d_number_groups = d_input->get<int>("number_groups");

  // Get relevant input parameters.
  if (input->check("outer_max_iters"))
    d_max_iters = input->get<int>("outer_max_iters");

  if (input->check("outer_tolerance"))
    d_tolerance = input->get<double>("outer_tolerance");

  if (input->check("outer_print_out"))
    d_print_out = input->get<int>("outer_print_out");

  if (input->check("outer_print_interval"))
    d_print_interval = input->get<int>("outer_print_interval");

  // We can turn off downscatter even if the material has
  // it and is set to use it.  This might be desirable when
  // we want to allow upscatter to be updated implicitly
  // in an outer eigenvalue iteration.
  if (material->downscatter()) d_downscatter = true;

  if (input->check("outer_downscatter"))
  {
    d_downscatter = input->get<int>("outer_downscatter");
  }

  // Create inner iteration.
  d_inner_solver = new SourceIteration<D>(input, state, mesh, material,
                                          quadrature, boundary, q_e, q_f);

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of GauseSeidel.cc
//---------------------------------------------------------------------------//

