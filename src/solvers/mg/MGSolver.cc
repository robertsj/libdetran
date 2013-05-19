//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolver.cc
 *  @brief MGSolver member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGSolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
MGSolver<D>::MGSolver(SP_state                  state,
                      SP_material               material,
                      SP_boundary               boundary,
                      const vec_externalsource &q_e,
                      SP_fissionsource          q_f,
                      bool                      multiply)
  : Base(state, material, boundary, q_e, q_f)
  , d_downscatter(false)
  , d_multiply(multiply)
{

  // Get relevant input parameters.
  if (d_input->check("outer_max_iters"))
    d_maximum_iterations = d_input->template get<int>("outer_max_iters");
  if (d_input->check("outer_tolerance"))
    d_tolerance = d_input->template get<double>("outer_tolerance");
  if (d_input->check("outer_print_level"))
    d_print_level = d_input->template get<int>("outer_print_level");
  if (d_input->check("outer_print_interval"))
    d_print_interval = d_input->template get<int>("outer_print_interval");
  if (d_input->check("adjoint"))
    d_adjoint = 0 != d_input->template get<int>("adjoint");

  // We can turn off downscatter even if the material has
  // it and is set to use it.  This might be desirable when
  // we want to allow upscatter to be updated implicitly
  // in an outer eigenvalue iteration.
  if (d_material->downscatter()) d_downscatter = true;
  if (d_input->check("outer_downscatter"))
  {
    d_downscatter = d_input->template get<int>("outer_downscatter");
  }

}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class MGSolver<_1D>;
template class MGSolver<_2D>;
template class MGSolver<_3D>;

} // end namespace detran



