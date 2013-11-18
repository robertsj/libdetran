//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenGD.cc
 *  @brief EigenGD member definitions.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenGD.hh"
#include "solvers/mg/MGSolverGMRES.hh"
#include "solvers/mg/MGDSA.hh"
#include "solvers/mg/MGCMDSA.hh"
#include "solvers/mg/MGTCDSA.hh"
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
  // By default, GD find largest; we want smallest, equal to
  // 1/k for largest k
  db->template put<int>("eigen_solver_which_value", 1);
  d_eigensolver = callow::EigenSolverCreator::Create(db);
  Assert(d_eigensolver);

  // Get callow solver parameter database
  SP_input pcdb;
  if (d_input->check("eigen_solver_pc_db"))
    pcdb = d_input->template get<SP_input>("eigen_solver_pc_db");
  else
    pcdb = NULL;
  //pcdb->display();
  d_eigensolver->set_operators(d_F, d_A, pcdb);

  // Preconditioner -- diffusion pc's assume k = 1.0 for now so
  // we need a nice way to update on the fly.  Moreover, the
  // diffusion PC's are not set up for boundaries.  The quick
  // fix is to allow the full vectors be passed to the
  // operators but to insert them in a possibly smaller vector
  // of correct size, thereby leaving the boundary terms unaccelerated
  // until a better approach is implemented
  typename MGSolverGMRES<D>::SP_pc pc;
  std::string pc_type = "default";
  if (d_input->check("eigen_solver_pc_type"))
    pc_type = d_input->template get<std::string>("eigen_solver_pc_type");
  MGTransportOperator<D> &TO = *mgs->get_operator();
  if (pc_type == "mgdsa")
  {
    pc = new MGDSA(d_input, d_material, d_mesh,
                   TO.sweepsource()->get_scatter_source(),
                   d_fissionsource,
                   0, true, d_adjoint);

  }
  else if (pc_type == "mgcmdsa")
  {
    pc = new MGCMDSA(d_input, d_material, d_mesh,
                     TO.sweepsource()->get_scatter_source(),
                     d_fissionsource,
                     0, true, d_adjoint);
    pc->build(1.0);
  }
  else if (pc_type == "mgtcdsa")
  {
    typename MGSolverGMRES<D>::SP_pc
      base_pc(new MGCMDSA(d_input,
                          d_material,
                          d_mesh,
                          TO.sweepsource()->get_scatter_source(),
                          d_fissionsource,
                          0,
                          true,
                          d_adjoint));

    pc = new MGTCDSA<D>(d_input,
                        d_material,
                        d_mesh,
                        TO.sweepsource()->get_scatter_source(),
                        d_fissionsource,
                        base_pc,
                        mgs->get_operator(),
                        0,
                        true,
                        d_adjoint);
    pc->build(1.0);
  }
  else
  {
    // use the default, which is a crude solve for (A-(1/k)F)v = r
  }

  if (pc) d_eigensolver->set_preconditioner(pc);

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
  d_state->set_eigenvalue(d_eigensolver->eigenvalue());
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
