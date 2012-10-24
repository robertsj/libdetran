//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   WGSolver.i.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  WGSolver inline member definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef detran_WGSOLVER_T_HH
#define detran_WGSOLVER_T_HH

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
  if (d_input->check("inner_print_out"))
    d_print_level = d_input->template get<int>("inner_print_level");
  if (d_input->check("inner_print_interval"))
    d_print_interval = d_input->template get<int>("inner_print_interval");

  //-------------------------------------------------------------------------//
  // MOMENT-TO-DISCRETE
  //-------------------------------------------------------------------------//

  SP_MtoD MtoD;
  // \todo Currently limited to 1 moment (0th order)
  MtoD = new detran_angle::MomentToDiscrete(0);
  MtoD->build(d_quadrature);

  //-------------------------------------------------------------------------//
  // SWEEP SOURCE
  //-------------------------------------------------------------------------//

  // Build the sweep source.
  d_sweepsource =
    new SweepSource<D>(d_state, d_mesh, d_quadrature, d_material,
                       MtoD, d_multiply);

  if (q_f) d_sweepsource->set_fission_source(q_f);
  Assert(!multiply or (multiply and q_f));

  // \todo Allow for discrete sources
  for (int i = 0; i < q_e.size(); ++i)
  {
    if (q_e[i]) d_sweepsource->set_moment_source(q_e[i]);
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
// 3-D
template <class D>
inline bool WGSolver<D>::set_sweep(std::string equation)
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
  return false;
}


} // end namespace detran

#endif /* detran_WGSOLVER_T_HH */

//---------------------------------------------------------------------------//
//              end of WGSolver.i.hh
//---------------------------------------------------------------------------//
