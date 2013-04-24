//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EigenDiffusion.cc
 *  @brief  EigenDiffusion
 *  @author Jeremy Roberts
 *  @date   Oct 26, 2012
 */
//---------------------------------------------------------------------------//

#include "EigenDiffusion.hh"
#include "MGDiffusionSolver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
EigenDiffusion<D>::EigenDiffusion(SP_mg_solver mg_solver)
  : Base(mg_solver)
{
  // Extract diffusion solver pointer
  MGDiffusionSolver<D>* mg_diff =
      dynamic_cast<MGDiffusionSolver<D>*>(mg_solver->solver().bp());
  Insist(mg_diff,
    "EigenDiffusion requires a multigroup diffusion solver");

  // Get mesh and materials
  d_mesh = d_mg_solver->mesh();
  d_material = d_mg_solver->material();

  // Create operators and vectors
  size_t problem_size = d_mesh->number_cells() * d_material->number_groups();
  d_phi  = new callow::Vector(problem_size, 0.0);
  d_work = new callow::Vector(problem_size, 1.0);

  // Get the loss operator and build the fission operator
  d_M = mg_diff->lossoperator();
  d_F = new DiffusionGainOperator(d_input, d_material, d_mesh, d_adjoint);

  d_F->print_matlab("eigF.out");
  d_M->print_matlab("eigM.out");


  // Create callow eigensolver and set operators
  SP_input db;
  if (d_input->check("eigen_solver_db"))
    db = d_input->template get<SP_input>("eigen_solver_db");
  d_eigensolver = Creator_T::Create(db);
  d_eigensolver->set_operators(d_F, d_M, db);
}

//---------------------------------------------------------------------------//
template <class D>
void EigenDiffusion<D>::solve()
{

  // solve the system
  d_eigensolver->solve(d_phi, d_work);

  // set the eigenvalue
  d_state->set_eigenvalue(d_eigensolver->eigenvalue());

  // fill the state flux vector
  int k = 0;
  for (int g = 0; g < d_material->number_groups(); g++)
    for (int i = 0; i < d_mesh->number_cells(); i++, k++)
      d_state->phi(g)[i] = (*d_phi)[k];

  // update the fission source (for use in rates, etc.)
  d_fissionsource->update();

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class EigenDiffusion<_1D>;
template class EigenDiffusion<_2D>;
template class EigenDiffusion<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file EigenDiffusion.cc
//---------------------------------------------------------------------------//
