//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenGD.cc
 *  @brief EigenGD member definitions.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenGD.hh"
#include "solvers/mg/MGSolverGMRES.hh"
#include "callow/solver/EigenSolverCreator.hh"

#include <iostream>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
EigenGD<D>::EigenGD(SP_mg_solver mg_solver)
  : Base(mg_solver)
{
  MGSolverGMRES<D>* mgs =
    dynamic_cast<MGSolverGMRES<D>*>(&(*d_mg_solver->solver()));
  Insist(mgs, "EigenGD requires GMRES for the MG problem to get the operator.");

  // Transport operator
  d_A = mgs->get_operator();
  d_A->sweeper()->set_update_boundary(false);

  // Fission operator
  d_F = new LHS_Operator_T(d_mg_solver);

  // Create vectors
  d_x  = new callow::Vector(d_A->number_rows(), 1.0);
  d_x0 = new callow::Vector(d_A->number_rows(), 1.0);

  // Get callow eigensolver parameter database
  SP_input db;
  if (d_input->check("eigen_solver_db"))
  {
    db = d_input->template get<SP_input>("eigen_solver_db");
  }
  else
  {
    db = new detran_utilities::InputDB();
  }
  db->template put<std::string>("eigen_solver_type", "gd");
  d_eigensolver = callow::EigenSolverCreator::Create(db);
  Assert(d_eigensolver);


//  d_F->compute_explicit("FF.out");
//  d_A->compute_explicit("AA.out");

  // Get callow solver parameter database
  if (d_input->check("eigen_solver_pc_db"))
    db = d_input->template get<SP_input>("eigen_solver_pc_db");
  else
    db = NULL;
  d_eigensolver->set_operators(d_F, d_A, db);

  // Preconditioner



}


//----------------------------------------------------------------------------//
template <class D>
void EigenGD<D>::solve()
{
  // Solve the problem
  d_eigensolver->solve(d_x, d_x0);

  // Fill the state
  size_t k = 0;
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_A->moments_size(); ++i, ++k)
    {
      d_state->phi(g)[i] = (*d_x)[k];
    }
  }

  // Fill the
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class EigenGD<_1D>;
template class EigenGD<_2D>;
template class EigenGD<_3D>;

} // end namespace detran


//----------------------------------------------------------------------------//
//              end of EigenGD.cc
//----------------------------------------------------------------------------//
