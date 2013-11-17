//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverGS.cc
 *  @brief MGSolverGS member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGSolverGS.hh"
#include <iostream>
#include <string>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
MGSolverGS<D>::MGSolverGS(SP_state                  state,
                          SP_material               material,
                          SP_boundary               boundary,
                          const vec_externalsource &q_e,
                          SP_fissionsource          q_f,
                          bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
  , d_lower(0)
  , d_lower_upscatter(d_material->upscatter_cutoff(d_adjoint))
  , d_upper(d_material->number_groups())
  , d_iterate(false)
  , d_norm_type("Linf")
{
  if (d_input->check("outer_norm_type"))
    d_norm_type = d_input->template get<std::string>("outer_norm_type");

  if ((!d_downscatter && d_maximum_iterations > 0 && d_number_groups > 1)
      || d_multiply)
  {
    d_iterate = true;
  }

  // For adjoint problems, the iteration is done in reverse.  This isn't
  // strictly required, but what is "upscatter iteration" in the forward mode
  // becomes "downscatter iteration" in the adjoint mode, and so this just
  // saves iterations.
  if (d_adjoint)
  {
    d_lower = d_number_groups - 1;
    d_lower_upscatter = d_material->upscatter_cutoff(d_adjoint);
    d_upper = -1;
  }

  // For multiplying problems, we assume iterations are all groups
  if (d_multiply) d_lower_upscatter = d_lower;


  Ensure(d_norm_type == "Linf" || d_norm_type == "L1" || d_norm_type == "L2");
}

//----------------------------------------------------------------------------//
template <class D>
int MGSolverGS<D>::number_sweeps() const
{
  return d_wg_solver->get_sweeper()->number_sweeps();
}

template class MGSolverGS<_1D>;
template class MGSolverGS<_2D>;
template class MGSolverGS<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of MGSolverGS.cc
//----------------------------------------------------------------------------//

