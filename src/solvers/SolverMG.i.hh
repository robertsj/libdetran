//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SolverMG.hh
 * \author robertsj
 * \date   Jun 19, 2012
 * \brief  SolverMG inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef MULTIGROUPSOLVER_I_HH_
#define MULTIGROUPSOLVER_I_HH_

// Detran
#include "detran_config.hh"
#include "SourceIteration.hh"
#ifdef DETRAN_ENABLE_PETSC
#include "InnerGMRES.hh"
#endif

namespace detran
{

template <class D>
SolverMG<D>::SolverMG(SP_input          input,
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
  , d_fissionsource(q_f)
  , d_downscatter(false)
  , d_max_iters(100)
  , d_tolerance(1e-5)
  , d_print_out(2)
  , d_print_interval(10)
  , d_multiply(false)
{
  Require(d_input);
  Require(d_state);
  Require(d_mesh);
  Require(d_material);
  Require(d_quadrature);
  Require(d_boundary);

  // Get relevant input parameters.
  if (input->check("outer_max_iters"))
    d_max_iters = input->get<int>("outer_max_iters");

  if (input->check("outer_tolerance"))
    d_tolerance = input->get<double>("outer_tolerance");

  if (input->check("outer_print_out"))
    d_print_out = input->get<int>("outer_print_out");

  if (input->check("outer_print_interval"))
    d_print_interval = input->get<int>("outer_print_interval");

  // Check if we need to include fission
  if (input->check("problem_type"))
  {
    if (input->get<std::string>("problem_type") == "multiply")
    {
      Insist(d_fissionsource,
            "Fission source must be constructed for multiplying problem");
      d_multiply = true;
    }
  }

  Assert(d_input->check("number_groups"));
  d_number_groups = d_input->get<int>("number_groups");

  // We can turn off downscatter even if the material has
  // it and is set to use it.  This might be desirable when
  // we want to allow upscatter to be updated implicitly
  // in an outer eigenvalue iteration.
  if (d_material->downscatter()) d_downscatter = true;
  if (input->check("outer_downscatter"))
  {
    d_downscatter = input->get<int>("outer_downscatter");
  }

  // Get the inner solver type and create.
  std::string inner_solver = "SI";
  if (input->check("inner_solver"))
  {
    inner_solver = input->get<std::string>("inner_solver");
  }
  if (inner_solver == "SI")
  {
    d_inner_solver = new SourceIteration<D>(input, state, mesh, material,
                                            quadrature, boundary, q_e, q_f);
  }
  else if (inner_solver == "GMRES")
  {
#ifdef DETRAN_ENABLE_PETSC
    d_inner_solver = new InnerGMRES<D>(input, state, mesh, material,
                                       quadrature, boundary, q_e, q_f);
#else
    THROW("InnerGMRES is not available because PETSc is not enabled.");
#endif
  }
  else
  {
    THROW("Unsupported inner solver type selected: "+inner_solver);
  }

}

} // end namespace detran


#endif /* MULTIGROUPSOLVER_I_HH_ */
