//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Eigensolver.cc
 *  @brief Eigensolver member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Eigensolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
Eigensolver<D>::Eigensolver(SP_mg_solver mg_solver)
  : d_mg_solver(mg_solver)
{
  Require(mg_solver);
  d_input = d_mg_solver->input();
  d_state = d_mg_solver->state();
  d_material = d_mg_solver->material();
  d_mesh  = d_mg_solver->mesh();
  d_fissionsource = d_mg_solver->fissionsource();
  Require(d_input);
  Require(d_state);
  Require(d_fissionsource);
  Require(d_input->check("number_groups"));
  d_number_groups = d_input->template get<int>("number_groups");

  // Get relevant input parameters.
  if (d_input->check("eigen_max_iters"))
    d_maximum_iterations = d_input->template get<int>("eigen_max_iters");
  if (d_input->check("eigen_tolerance"))
    d_tolerance = d_input->template get<double>("eigen_tolerance");
  if (d_input->check("eigen_print_level"))
    d_print_level = d_input->template get<int>("eigen_print_level");
  if (d_input->check("eigen_print_interval"))
    d_print_interval = d_input->template get<int>("eigen_print_interval");
}

//----------------------------------------------------------------------------//
template <class D>
Eigensolver<D>::~Eigensolver()
{
  /* ... */
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATION
//----------------------------------------------------------------------------//

template class Eigensolver<_1D>;
template class Eigensolver<_2D>;
template class Eigensolver<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of Eigensolver.cc
//----------------------------------------------------------------------------//
