/*
 * TransientSolver.cc
 *
 *  Created on: Sep 6, 2020
 *      Author: rabab
 */

#include "TransientSolver.hh"


TransientSolver::TransientSolver(SP_input inp, SP_mesh mesh, SP_material material,
		                         SP_matrix flux_basis, SP_matrix precursors_basis)
:d_mesh(mesh),
 d_material(material),
 d_inp(inp),
 d_precursors_basis(precursors_basis),
 d_flux_basis(flux_basis)
{
 d_num_cells = d_mesh->number_cells();
 d_number_groups = d_material->number_groups();
 d_precursors_group = d_material->number_precursor_groups();
 d_rf = d_flux_basis->number_columns();
 d_rc = d_precursors_basis->number_columns();


 if (d_inp->check("ts_final_time"))
    d_final_time = d_inp->template get<double>("ts_final_time");
  Assert(d_final_time > 0.0);
  if (d_inp->check("ts_step_size"))
    d_dt = d_inp->template get<double>("ts_step_size");
  Assert(d_dt > 0.0);

  // Compute the number of steps.  May result in longer time than requested!
  d_number_steps = std::ceil(d_final_time / d_dt);

  d_P0 = new callow::Vector(d_precursors_group*d_num_cells, 0.0);
  d_phi0_r = new callow::Vector(d_rf, 0.0);
  d_P0_r = new callow::Vector (d_rc, 0.0);
  d_phi0 = new callow::Vector(d_num_cells* d_number_groups, 0.0);

  d_fluxes = new callow::MatrixDense(d_number_groups*d_num_cells, d_number_steps);
  d_precursors = new callow::MatrixDense(d_num_cells*d_precursors_group, d_number_steps);

  d_solver = LinearSolverCreator::Create(d_inp);

}

//------------------------------------------------------------------------------------//

void TransientSolver::initialize_precursors()
{
 d_fissionsource = new FissionSource(d_state, d_mesh, d_material);
 d_fissionsource->update();

 const State::moments_type &fd = d_fissionsource->density();

 const vec_int &mt = d_mesh->mesh_map("MATERIAL");

 for (int i = 0; i < d_precursors_group; ++i)
 {
  double inv_lambda = 1.0 / d_material->lambda(i);
  for (int cell = 0; cell < d_num_cells; ++cell)
 {
  (*d_P0)[cell + i*d_num_cells] = inv_lambda * d_material->beta(mt[cell], i) * fd[cell];
 }

 }

}
//------------------------------------------------------------------------------------//

void TransientSolver::ProjectInitial()
{
 for (int g=0; g<d_number_groups; g++)
 {
  for (int cell=0; cell<d_num_cells; cell++)
  {
   (*d_phi0)[cell + g*d_num_cells] = d_state->phi(g)[cell];
  }
 }

 d_flux_basis->multiply_transpose(*d_phi0, *d_phi0_r);

 //project initial precursors
 d_precursors_basis->multiply_transpose(*d_P0, *d_P0_r);
}

//------------------------------------------------------------------------------------//

void TransientSolver::UpdatePrecursors(double dt)
{
 //update reduced precursors
 callow::Vector d_P_0(d_rc, 0.0);

 d_F->multiply(*d_phi_r, d_P_0);

 d_P_0.scale(dt);

 d_P_0.add(d_P0_r);
 // scale dP_0 by dt
 LinearSolverCreator::SP_solver solver;
 // solve Ax = b  , need callow solver
 solver = LinearSolverCreator::Create(d_inp);

 solver->set_operators(d_P);

 d_P_r = new callow::Vector(d_rc, 0.0);
 solver->solve(d_P_0, *d_P_r);
 std::cout <<" pp "<< (*d_P_r)[0] << "\n";
 std::cout << "pp " << (*d_P_r)[1] << "\n";

}

//------------------------------------------------------------------------------------//

void TransientSolver::Solve(SP_state initial_state)
{
 bool rom = true;
 //d_material->update(0.0, 0, 1, false);

 d_state = initial_state;

 initialize_precursors();

 ProjectInitial();

 double t = 0.0;

 for (int i=0 ; i< d_number_steps; i++)
 {
  std::cout << "step " << i << "\n";
  t += d_dt;
  step(t, d_dt, i+1,  rom);
 }

 d_precursors->print_matlab("d_precursors.txt");
 d_fluxes->print_matlab("d_flux.txt");

 std::cout <<  "################### END 1 ############################ \n";
}

//------------------------------------------------------------------------------------//

void TransientSolver::UpdateOperator(double dt)
{
 // get the matrices
 Kinetic_Mat K(d_inp, d_mesh, d_material, d_flux_basis, d_precursors_basis);

 // need to update the material here
 int rf = d_flux_basis->number_columns();
 int rc = d_precursors_basis->number_columns();

 d_P = new callow::MatrixDense(rc, rc);
 d_P = K.Mat1();

 d_D = new callow::MatrixDense(rf, rc);
 d_D = K.Mat2();

 d_F = new callow::MatrixDense(rc, rf);
 d_F = K.Mat3();

 // if diffusion
 SP_lossoperator L (new DiffusionLossOperator(d_inp, d_material, d_mesh, false, 0.0, false, 1.0));
 std::cout << "sigma "<< d_material->sigma_t(2,1) << "\n";
 std::cout << "chi "<< d_material->chi(1, 0) << "\n";
 std::cout << "chi "<< d_material->chi(1, 1) << "\n";



 SP_gainoperator G (new DiffusionGainOperator(d_inp, d_material, d_mesh, false));

 OperatorProjection Projector(1);

 d_L = new callow::MatrixDense(rf, rf);

 Projector.SetOperators(L, d_flux_basis);
 Projector.Project(d_L);
 //d_L->print_matlab("d_L.txt");
 printf(" %1.16f \n", (*L)(0, 0));
 printf(" %1.16f \n", (*L)(1, 1));
 printf(" %1.16f \n", (*L)(2, 2));
 printf(" %1.16f \n", (*L)(2, 0));

 //std::cout << "L " << (*d_L)[4, 5] << "\n";
 //std::cout << "L " << (*d_L)[6, 9] << "\n";
 //d_L->print_matlab("d_L");

 d_G = new callow::MatrixDense(rf, rf);
 Projector.SetOperators(G, d_flux_basis);
 Projector.Project(d_G);

 //compute (I + dt*Pr)
  for (int i=0; i<d_rc; i++)
  {
   for (int j=0; j<d_rc; j++)
   {
 	(*d_P)(i, j) = (*d_P)(i, j)*dt;
   }
   // add the identity matrix
   d_P->insert(i, i, 1.0, 1);
  }

  d_solver->set_operators(d_P);

  // need to compute P^-1 * F -> PF
  MatrixDense PF(d_rc, d_rf);
  for (int i=0; i < d_rf; i++)
  {
   Vector y(d_rc, 0.0);
   Vector x(d_rc, 0.0);
   for (int j=0; j< d_rc; j++)
   {
    x[j] = (*d_F)(j, i);
   }
   d_solver->solve(x, y);
   for (int k=0; k<d_rc; k++)
   {
    PF.insert(k, i, y[k]);
   }
  }

  PF.print_matlab("PF.txt");

  /*
  std::cout << "P 0 0  " << (*d_P)(0, 0);
  std::cout << "P 1 0  " << (*d_P)(1, 0);
  std::cout << "P 3 2  " << (*d_P)(3, 2);
  std::cout << "P 4 4  " << (*d_P)(4, 4);

  std::cout << "*******************************\n";
  std::cout << "F 4 7 " << (*d_F)(4,7) << "\n";
  std::cout << "F 4 3 " << (*d_F)(4,3) << "\n";
  std::cout << "PF 0 0  " << PF(0,0) << "\n";
  std::cout << "PF 1 0  " << PF(1,0) << "\n";
  std::cout << "PF 2 3  " << PF(2,3) << "\n";
  std::cout << "PF 3 4  " << PF(3,4) << "\n";
  std::cout << "PF 4 4  " << PF(4,4) << "\n";
  std::cout << "PF 4 7  " << PF(4,7) << "\n";
  std::cout << "PF 4 8  " << PF(4,8) << "\n";
  */

 // compute D@PF =  B
 d_B = new callow::MatrixDense(rf, rf);
 for (int i=0; i<d_rf ;i++)
 {
  for (int j=0; j<d_rf; j++)
  {
   for (int k=0; k<d_rc; k++)
   {
    double value = (*d_D)(i, k)*PF(k,j);
    d_B->insert(i, j, value, 1);
   }
  }
 }

 d_B->print_matlab("B.txt");

 // std::cout << "DPF 0 0  " << (*d_B)(0,0) << "\n";
 // std::cout << "DPF 1 0  " << (*d_B)(1,0) << "\n";
 // std::cout << "DPF 3 4  " << (*d_B)(3,4) << "\n";
 // std::cout << "DPF 5 7  " << (*d_B)(4,7) << "\n";

 d_A = new callow::MatrixDense(rf, rf);
 // compute the matrix
 for (int i=0; i<d_rf; i++)
 {
  for (int j=0; j<d_rf; j++)
  {
   double value = (*d_L)(i, j) - (*d_G)(i, j) - (*d_B)(i, j)*dt;
   d_A->insert(i, j, value, 1);
  }
 }

 // std::cout << "A 0 0 " << (*d_A)(0, 0) << "\n";
 // std::cout << "A 1 0 " << (*d_A)(1, 0) << "\n";
 // std::cout << "A 1 1 " << (*d_A)(1, 1) << "\n";
 // std::cout << "A 3 4 " << (*d_A)(3, 4) << "\n";
}
//------------------------------------------------------------------------------------//


void TransientSolver::compute_b(double d_dt)
{

 d_solver->set_operators(d_P);
 // compute P^-1*c0_r

 Vector Pc0(d_rc, 0.0);

 // need to change d_P0
 d_solver->solve(*d_P0_r, Pc0);

 // compute DP^-1*c0
 Vector DPc(d_rf, 0.0);

//DPc =  new Vector(rf, 0.0);
 d_D->multiply(Pc0, DPc);

 // compute b
 d_b = new callow::Vector (d_rf, 0.0);

 int rfg = d_rf/d_number_groups;

 for (int i=0; i<rfg; i++)
 {
  for(int g=0; g<d_number_groups; g++)
  {
  (*d_b)[i + g*rfg] =
		  DPc[i + g*rfg] + ((*d_phi0_r)[i + g*rfg]/(d_material->velocity(g) * d_dt));
  //std::cout << "b " << (*d_b)[i + g*rfg] << "\n";
  }
 }
}

void TransientSolver::step(double t, double dt, int step_idx, bool rom_flag)
{
 // need to update the material here
 d_material->update(t, dt, 1, true);

 UpdateOperator(dt);

 compute_b(dt);

 d_solver->set_operators(d_A);
 d_phi_r = new callow::Vector(d_rf, 0.0);
 d_solver->solve(*d_b, *d_phi_r);

 UpdatePrecursors(dt);

 // project back and store in d_phi if the final step

 d_flux_basis->multiply(*d_phi_r, *d_phi0);

 double phi[d_num_cells*d_number_groups];

 for (int i=0; i < d_num_cells*d_number_groups; i++)
 {
  phi[i] = (*d_phi0)[i];
 }

d_fluxes->insert_col(step_idx, phi, 0);

 d_precursors_basis->multiply(*d_P_r, *d_P0);
 double precursors[d_num_cells*d_precursors_group];



 for (int i =0; i <d_num_cells*d_precursors_group; i++)
 {
  precursors[i] = (*d_P0)[i] ;
 }
 d_precursors->insert_col(step_idx,  precursors);

 // cycle
 *d_phi0_r = *d_phi_r;

 *d_P0_r = *d_P_r;

}


