//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionEigensolver.cc
 * \author robertsj
 * \date   Jul 26, 2012
 * \brief  DiffusionEigensolver class definition.
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.hh"

// Diffusion
#include "DiffusionEigensolver.hh"

// System
#include <string>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <typeinfo>

namespace detran
{

template <class D>
DiffusionEigensolver<D>::DiffusionEigensolver(SP_input    input,
                                              SP_material material,
                                              SP_mesh     mesh,
                                              SP_state    state)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_state(state)
  , d_boundary(new BoundaryDiffusion<D>(input, mesh))
  , d_adjoint(false)
  , d_solver_type("power")
  , d_maximum_iterations(10000)
  , d_tolerance(1e-6)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_state);
  Require(d_boundary);

  d_problem_size = d_mesh->number_cells() * d_material->number_groups();

  // Select solver type
  if (d_input->check("diffusion_eigensolver_type"))
  {
    d_solver_type =
        d_input->template get<std::string>("diffusion_eigensolver_type");
  }
  // Select solve adjoint
  if (d_input->check("diffusion_adjoint"))
  {
    d_adjoint = d_input->template get<int>("diffusion_adjoint");
  }
  // maximum iterations
  if (d_input->check("diffusion_maximum_iterations"))
  {
    d_maximum_iterations =
      d_input->template get<int>("diffusion_maximum_iterations");
  }
  // tolerance
  if (d_input->check("diffusion_tolerance"))
  {
    d_tolerance =
      d_input->template get<double>("diffusion_tolerance");
  }

  // Create operators and vectors
  d_phi  = new Vector_T(d_problem_size, 0.0);
  d_work = new Vector_T(d_problem_size, 1.0);
  d_F = new DiffusionGainOperator(d_input, d_material, d_mesh, d_adjoint);
  d_M = new DiffusionLossOperator(d_input, d_material, d_mesh,
                                  false, d_adjoint);
  d_F->print_matlab("eigF.out");
  d_M->print_matlab("eigM.out");

  // Create solver
  d_solver = Creator_T::Create(d_solver_type,
                               d_tolerance,
                               d_maximum_iterations);

  // Set the operators
  d_solver->set_operators(d_F, d_M);
}

template <class D>
void DiffusionEigensolver<D>::DiffusionEigensolver::solve()
{

  // Solve the system
  d_solver->solve(d_phi, d_work);

  // set the eigenvalue
  d_state->set_eigenvalue(d_solver->eigenvalue());
  // fill the state flux vector
  int k = 0;
  for (int g = 0; g < d_material->number_groups(); g++)
    for (int i = 0; i < d_mesh->number_cells(); i++, k++)
      d_state->phi(g)[i] = (*d_phi)[k];

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class DiffusionEigensolver<_1D>;
template class DiffusionEigensolver<_2D>;
template class DiffusionEigensolver<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file DiffusionEigensolver.cc
//---------------------------------------------------------------------------//
