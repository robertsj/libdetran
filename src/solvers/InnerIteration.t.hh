//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIteration.i.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerIteration inline member definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef INNERITERATION_T_HH_
#define INNERITERATION_T_HH_

#include <iostream>

namespace detran
{

template <class D>
InnerIteration<D>::InnerIteration(SP_input          input,
                                  SP_state          state,
                                  SP_mesh           mesh,
                                  SP_material       material,
                                  SP_quadrature     quadrature,
                                  SP_boundary       boundary,
                                  SP_externalsource q_e,
                                  SP_fissionsource  q_f)
  :  d_input(input)
  ,  d_state(state)
  ,  d_mesh(mesh)
  ,  d_material(material)
  ,  d_quadrature(quadrature)
  ,  d_boundary(boundary)
  ,  d_max_iters(100)
  ,  d_tolerance(1e-5)
  ,  d_print_out(2)
  ,  d_print_interval(10)
{
  Require(input);
  Require(state);
  Require(mesh);
  Require(material);
  Require(quadrature);
  Require(boundary);

  // Get relevant input parameters.
  if (input->check("inner_max_iters"))
  {
    d_max_iters = input->get<int>("inner_max_iters");
  }
  if (input->check("inner_tolerance"))
  {
    d_tolerance = input->get<double>("inner_tolerance");
  }
  if (input->check("inner_print_out"))
  {
    d_print_out = input->get<int>("inner_print_out");
  }
  if (input->check("inner_print_interval"))
  {
    d_print_interval = input->get<int>("inner_print_interval");
  }

  // Moments-to-Discrete
  SP_MtoD MtoD;
  MtoD = new MomentToDiscrete<D>(0); // 1 moment (0th order)
  MtoD->build(quadrature);

  // Build the sweep source.
  d_sweepsource =
      new SweepSource<D>(state, mesh, quadrature, material, MtoD);

  if (q_f) d_sweepsource->set_fission_source(q_f);

  if (q_e)
  {
    // \todo Come up with a nice way to add/remove sources.  Perhaps
    //       also allow names?
    d_sweepsource->set_moment_source(q_e);
  }

  // Get equation type.  Diamond-difference is default.
  std::string equation = "dd";
  if (d_input->check("equation"))
  {
    equation = d_input->get<std::string>("equation");
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

// 3-D
template <class D>
inline bool InnerIteration<D>::set_sweep(std::string equation)
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

// 2-D
template <>
inline bool InnerIteration<_2D>::set_sweep(std::string equation)
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

// 1-D
template <>
inline bool InnerIteration<_1D>::set_sweep(std::string       equation)
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

#endif /* INNERITERATION_T_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIteration.i.hh
//---------------------------------------------------------------------------//
