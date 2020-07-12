/*
 * ROM_Manager.cc
 *
 *  Created on: Jul 8, 2020
 *      Author: rabab
 */

#include "ROM_Manager.hh"

template <class D>
ROM_Manager<D>::ROM_Manager(SP_input inp, SP_mesh mesh, SP_material mat, std::string operator_type)
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

template <class D>
void ROM_Manager<D>::Set_FullOperators()
{
	if (d_operator_type == "diffusion")

	{
	  d_A = new DiffusionLossOperator(d_input, d_mat, d_mesh, false, 0.0, false, 1.0);
	 // d_A->print_matlab("A.txt");
	  d_B = new DiffusionGainOperator(d_input, d_mat, d_mesh, false);
	  //d_B->print_matlab("B.txt");
	}

	else if (d_operator_type == "EnergyDependent")
	{
	  d_input->put<std::string>("outer_solver", "GMRES");
	  d_input->put<std::string>("eigen_solver", "GD");
	  SP_mg_solver mg_solver_ED;
	  mg_solver_ED = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, true, true);
	  mg_solver_ED->setup();
	  mg_solver_ED->set_solver();
	  d_A = new LHS_Operator_T(mg_solver_ED);
	  //d_A->compute_explicit("EnergyDependent");
	}

	else if (d_operator_type == "EnergyIndependent")
	{
	 SP_mg_solver mg_solver;
	 mg_solver = new FixedSourceManager<_1D>(d_input, d_mat, d_mesh, false, true);
	 mg_solver->setup();
	 mg_solver->set_solver();
	 d_A = new Operator_T(mg_solver);
     std::cout << d_A->number_columns() << "\n";
	 //d_A->compute_explicit("/home/rabab/Desktop/EnergyIndependent");
	}

}

template <class D>
void ROM_Manager<D>::Solve(SP_matrix d_U)
{
	detran::ProjectedOperator P(1);

	P.SetOperators(d_A, d_U);

    int r = d_U->number_columns();

	SP_matrix Ar;
	Ar = new callow::MatrixDense(r, r);
	P.Project(Ar);
	//Ar->print_matlab("Ar.txt");

	SP_eigensolver eigensolver;
	eigensolver = Creator_T::Create(d_input);

	if (d_operator_type == "diffusion")
	{
		P.SetOperators(d_B, d_U);
		SP_matrix Br;
		Br = new callow::MatrixDense(r, r);
		P.Project(Br);
		//Br->print_matlab("Br.txt");
		eigensolver->set_operators(Br, Ar);
	}

	else eigensolver->set_operators(Ar);

    SP_vector x_rom;
    SP_vector x1;
    x_rom = new callow::Vector(r, 0.0);
    x1 = new callow::Vector(r, 1.0);
    eigensolver->solve(x_rom, x1);

    std::cout << (*x_rom)[0] << "***********\n";

    d_keff =  eigensolver->eigenvalue();

    // reconstruct
    int n = d_U->number_rows();
    SP_vector x_fom;
    x_fom = new callow::Vector(n, 0.0);
    d_U->multiply(x_rom, x_fom);
    }

SOLVERS_INSTANTIATE_EXPORT(ROM_Manager<_1D>)
SOLVERS_INSTANTIATE_EXPORT(ROM_Manager<_2D>)
SOLVERS_INSTANTIATE_EXPORT(ROM_Manager<_3D>)
