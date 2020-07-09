/*
 * ROM_Manager.cc
 *
 *  Created on: Jul 8, 2020
 *      Author: rabab
 */

#include "ROM_Manager.hh"

ROM_Manager::ROM_Manager(SP_input inp, SP_mesh mesh, SP_material mat, std::string operator_type)
:d_input(inp),
 d_mesh(mesh),
 d_mat(mat),
 d_operator_type(operator_type)

{
 Require(inp);
 Require(mat);
 Require(mesh);
 ROM_Manager::Set_FullOperators();
}


void ROM_Manager::SetBasis()
{

}


void ROM_Manager::Set_FullOperators()
{
	if (d_operator_type ==  "diffusion")
	{
		SP_lossoperator d_A (new detran::DiffusionLossOperator(d_input, d_mat, d_mesh, false, 0.0, false, 1.0));
		SP_gainoperator d_B (new detran::DiffusionGainOperator(d_input, d_mat, d_mesh, false));
	}

	else if (d_operator_type == "EnergyDependent")
	{
		  d_input->put<std::string>("outer_solver", "GMRES");
		  d_input->put<std::string>("eigen_solver", "GD");
		  SP_mg_solver mg_solver_ED;
		  mg_solver_ED = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
		  mg_solver_ED->setup();
		  mg_solver_ED->set_solver();

		   d_A = new LHS_Operator_T(mg_solver_ED);
	}

	else if (d_operator_type == "EnergyIndependent")
	{
		d_input->put<std::string>("outer_solver", "GMRES");
		d_input->put<std::string>("eigen_solver", "GD");
		SP_mg_solver mg_solver;
		mg_solver = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
		mg_solver->setup();
	    mg_solver->set_solver();
	    d_A = new Operator_T(mg_solver);
	}

}



void ROM_Manager::Solve(SP_matrix d_U)
{


	detran::ProjectedOperator P(1);
	P.SetOperators(d_A, d_U);

    int r = d_U->number_rows();

	SP_matrix Ar;

	P.Project(Ar);
	Ar = new callow::MatrixDense(r, r);

    SP_eigensolver eigensolver;
    eigensolver = Creator_T::Create(d_input);
    eigensolver->set_operators(Ar);

    SP_vector x_rom;
    SP_vector x1;
    x_rom = new callow::Vector(r, 0.0);
    x1 = new callow::Vector(r, 1.0);
    eigensolver->solve(x_rom, x1);

    keff =  eigensolver->eigenvalue();
    int n = d_U->number_columns();

    // reconstruct
    SP_vector x_fom;
    x_fom = new callow::Vector(n, 0.0);
    d_U->multiply(x_rom, x_fom);

}
