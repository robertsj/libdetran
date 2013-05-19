//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  WGSolver.cc
 *  @brief WGSolver member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "WGSolver.hh"
#include "transport/Sweeper1D.cc"
#include "transport/Sweeper2D.cc"
#include "transport/Sweeper3D.cc"
#include "transport/Sweeper2DMOC.cc"
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
WGSolver<D>::WGSolver(SP_state                  state,
                      SP_material               material,
                      SP_quadrature             quadrature,
                      SP_boundary               boundary,
                      const vec_externalsource &q_e,
                      SP_fissionsource          q_f,
                      bool                      multiply)
  : Base(state, material, boundary, q_e, q_f)
  , d_quadrature(quadrature)
  , d_multiply(multiply)
{
  // Preconditions
  Require(quadrature);

  //-------------------------------------------------------------------------//
  // SET TOLERANCES
  //-------------------------------------------------------------------------//

  if (d_input->check("inner_max_iters"))
    d_maximum_iterations = d_input->template get<int>("inner_max_iters");
  if (d_input->check("inner_tolerance"))
    d_tolerance = d_input->template get<double>("inner_tolerance");
  if (d_input->check("inner_print_level"))
    d_print_level = d_input->template get<int>("inner_print_level");
  if (d_input->check("inner_print_interval"))
    d_print_interval = d_input->template get<int>("inner_print_interval");

  //-------------------------------------------------------------------------//
  // MOMENT-TO-DISCRETE
  //-------------------------------------------------------------------------//

  SP_MtoD MtoD;
  // \todo Currently limited to 1 moment (0th order)
  MtoD = new detran_angle::MomentToDiscrete(d_state->get_momentindexer());
  MtoD->build(d_quadrature);

  //-------------------------------------------------------------------------//
  // SWEEP SOURCE
  //-------------------------------------------------------------------------//

  // Build the sweep source.
  d_sweepsource =
    new SweepSource<D>(d_state, d_mesh, d_quadrature, d_material,
                       MtoD, d_multiply);

  if (q_f) d_sweepsource->set_fission_source(q_f);
  Assert(!multiply || (multiply && q_f));

  // Add any moment and discrete sources present.
  for (size_t i = 0; i < q_e.size(); ++i)
  {
    if (q_e[i])
    {
      if (!q_e[i]->is_discrete())
        d_sweepsource->set_moment_source(q_e[i]);
      else
        d_sweepsource->set_discrete_source(q_e[i]);
    }
  }

  //-------------------------------------------------------------------------//
  // SWEEPER
  //-------------------------------------------------------------------------//

  // Get equation type.  Diamond-difference is default.
  std::string equation = "dd";
  if (d_input->check("equation"))
  {
    equation = d_input->template get<std::string>("equation");
  }

  // Set the equation.
  if (!set_sweep(equation))
  {
    std::string message = "Equation type ";
    message.append(equation);
    message.append(" is not supported for this dimension.");
    THROW(message);
  }
}

//---------------------------------------------------------------------------//
// Default
template <class D>
inline bool WGSolver<D>::set_sweep(std::string equation)
{
  THROW("NOT IMPLEMENTED");
}

//---------------------------------------------------------------------------//
// 3-D
template <>
inline bool WGSolver<_3D>::set_sweep(std::string equation)
{
  if (equation == "dd")
  {
    d_sweeper = new Sweeper3D<Equation_DD_3D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  return false;
}

//---------------------------------------------------------------------------//
// 2-D
template <>
inline bool WGSolver<_2D>::set_sweep(std::string equation)
{

  if (equation == "dd")
  {
    d_sweeper = new Sweeper2D<Equation_DD_2D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  else if (equation == "sc")
  {
    d_sweeper = new Sweeper2D<Equation_SC_2D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  else if (equation == "sd")
  {
    d_sweeper = new Sweeper2D<Equation_SD_2D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  else if (equation == "scmoc")
  {
    d_sweeper = new Sweeper2DMOC<Equation_SC_MOC>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }

  return false;
}

//---------------------------------------------------------------------------//
// 1-D
template <>
inline bool WGSolver<_1D>::set_sweep(std::string       equation)
{
  if (equation == "dd")
  {
    d_sweeper = new Sweeper1D<Equation_DD_1D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  else if (equation == "sd")
  {
    d_sweeper = new Sweeper1D<Equation_SD_1D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  else if (equation == "sc")
  {
    d_sweeper = new Sweeper1D<Equation_SC_1D>(
      d_input, d_mesh, d_material, d_quadrature, d_state, d_boundary,
      d_sweepsource);
    return true;
  }
  return false;
}


//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class WGSolver<_1D>;
template class WGSolver<_2D>;
template class WGSolver<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of WGSolver.cc
//---------------------------------------------------------------------------//
