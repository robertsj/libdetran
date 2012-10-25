//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Eigensolver.cc
 *  @author robertsj
 *  @date   Jun 18, 2012
 *  @brief  Eigensolver class definition.
 */
//---------------------------------------------------------------------------//

#include "Eigensolver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
Eigensolver<D>::Eigensolver(SP_mg_solver mg_solver)
  : d_mg_solver(mg_solver)
{
  // Preconditions
  Require(mg_solver);
  d_input = d_mg_solver->input();
  d_state = d_mg_solver->state();
  d_fissionsource = d_mg_solver->fissionsource();
  Require(d_input);
  Require(d_state);
  Require(d_fissionsource);

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

//---------------------------------------------------------------------------//
template <class D>
Eigensolver<D>::~Eigensolver()
{
  /* ... */
}

//---------------------------------------------------------------------------//
// Explicit instantiations
template class Eigensolver<_1D>;
template class Eigensolver<_2D>;
template class Eigensolver<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Eigensolver.cc
//---------------------------------------------------------------------------//



