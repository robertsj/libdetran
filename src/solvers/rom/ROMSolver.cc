/*
 * ROM_Manager.cc
 *
 *  Created on: Jul 8, 2020
 *      Author: rabab
 */

#include "ROMSolver.hh"

template <class D>
ROMSolver<D>::ROMSolver(SP_input inp, SP_mesh mesh, SP_material mat)
:d_input(inp),
 d_mesh(mesh),
 d_mat(mat)
{
 Require(inp);
 Require(mat);
 Require(mesh);

 // not sure what input can tell that we want energyindependent operator,
 // so will let it the default until Roberts give his opinion
 d_operator = "EnergyDependent";

 if (d_input->check("equation"))
 {
   if (d_input->get<std::string>("equation") == "diffusion")
    d_operator = "diffusion";

    else if (d_input->get<std::string>("equation") == "dd")
	  d_operator = "EnergyIndependent";
 }

 ROMSolver::Set_FullOperators();
}

template <class D>
void ROMSolver<D>::Set_FullOperators()
{
  if (d_operator == "diffusion")
  {
    SP_lossoperator A (new DiffusionLossOperator(d_input, d_mat, d_mesh, false, 0.0, false, 1.0));
	SP_gainoperator B (new DiffusionGainOperator(d_input, d_mat, d_mesh, false));

	d_A = A;
	d_B = B;
  }

  else if (d_operator == "EnergyDependent")
  {
    d_input->put<std::string>("outer_solver", "GMRES");
	d_input->put<int>("outer_krylov_group_cutoff", 0);
	d_input->put<std::string>("eigen_solver", "GD");
	SP_mg_solver mg_solver;
	mg_solver = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
	mg_solver->setup();
	mg_solver->set_solver();
	d_B = new LHS_Operator_T(mg_solver);

	typename RHS_Operator_T::SP_operator A;
	MGSolverGMRES<D>* mgs =
	dynamic_cast<MGSolverGMRES<D>*>(&(*mg_solver->solver()));
	Insist(mgs, "EigenGD requires GMRES for the MG problem to get the operator.");

	// Transport operator
	A = mgs->get_operator();
	A->sweeper()->set_update_boundary(false);
	d_A = A;
  }

  else if (d_operator == "EnergyIndependent")
  {
   SP_mg_solver mg_solver;
   mg_solver = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
   mg_solver->setup();
   mg_solver->set_solver();
   d_A = new Operator_T(mg_solver);
  }
}

template <class D>
void ROMSolver<D>::Solve(SP_matrix d_U, SP_vector sol)
{
  detran::OperatorProjection P(1);

  P.SetOperators(d_A, d_U);

  int d_r = d_U->number_columns();

  SP_matrix Ar;
  Ar = new callow::MatrixDense(d_r, d_r);
  P.Project(Ar);

  SP_eigensolver eigensolver;
  eigensolver = Creator_T::Create(d_input);

  if (d_operator == "diffusion" || d_operator == "EnergyDependent")
  {
   P.SetOperators(d_B, d_U);
   SP_matrix Br;
   Br = new callow::MatrixDense(d_r, d_r);
   P.Project(Br);
   eigensolver->set_operators(Br, Ar);
  }

  else eigensolver->set_operators(Ar);

  SP_vector x_rom;
  SP_vector x;
  x_rom = new callow::Vector(d_r, 0.0);
  x = new callow::Vector(d_r, 1.0);
  eigensolver->solve(x_rom, x);

  d_keff =  eigensolver->eigenvalue();

  // reconstruct
  int d_n = d_U->number_rows();

  d_U->multiply(x_rom, sol);

  // correct the direction of the vector if negative
  for (int i=0; i<d_n; i++)
  {
    if (signbit((*sol)[i]))
	{
      (*sol)[i] *= -1;
	}
  }

  }


//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

SOLVERS_INSTANTIATE_EXPORT(ROMSolver<_1D>)
SOLVERS_INSTANTIATE_EXPORT(ROMSolver<_2D>)
SOLVERS_INSTANTIATE_EXPORT(ROMSolver<_3D>)
