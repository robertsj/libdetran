/*
 * TransientSolver.cc
 *
 *  Created on: Sep 6, 2020
 *      Author: rabab
 */

#include "TransientSolver.hh"


TransientSolver::TransientSolver(SP_input inp, SP_mesh mesh, SP_material material, SP_matrix flux_basis, SP_matrix precursors_basis)
:d_mesh(mesh),
 d_material(material),
 d_inp(inp),
 d_precursors_basis(precursors_basis),
 d_flux_basis(flux_basis),
 d_dt(1)
{
 d_num_cells = d_mesh->number_cells();
 d_number_groups = d_material->number_groups();
 d_precursors_group = d_material->number_precursor_groups();
 d_rf = d_flux_basis->number_columns();
 d_rc = d_precursors_basis->number_columns();
 d_fissionsource = new FissionSource(d_state, d_mesh, d_material);
}


void TransientSolver::initialize_precursors()
{
  d_fissionsource->update();
  const State::moments_type &fd = d_fissionsource->density();
  const vec_int &mt = d_mesh->mesh_map("MATERIAL");
  for (int i = 0; i < d_material->number_precursor_groups(); ++i)
  {
   double inv_lambda = 1.0 / d_material->lambda(i);
   for (int cell = 0; d_num_cells; ++cell)
  {
   d_precursor->C(i)[cell] =
   inv_lambda * d_material->beta(mt[cell], i) * fd[cell];
   //printf("%16.9f %16.9f %16.9f %16.9f \n", fd[cell], inv_lambda, d_material->beta(mt[cell], i), d_precursor->C(i)[cell]);
  }
  }
 // multiply_transpose function needs the input vector to be callow vector type,
 // want to check if there is a better way.
  int p_size = d_precursors_group*d_num_cells;

  d_P0 = new callow::Vector(p_size, 0.0);
  for (int i = 0; i < d_material->number_precursor_groups(); ++i)
  {
   for (int cell = 0; cell < d_num_cells; ++cell)
   {
     (*d_P0)[cell + i*d_num_cells] = d_precursor->C(i)[cell];
   }
  }

}


void TransientSolver::ProjectInitial(SP_state initial_state)
{
 int flux_size = d_num_cells* d_number_groups;
 d_phi0 = new callow::Vector(flux_size, 0.0);

 for (int g=0; g<d_number_groups; g++)
 {
  for (int cell=0; cell<d_num_cells; cell++)
  {
   (*d_phi0)[cell + g*d_num_cells] = initial_state->phi(g)[cell];
  }
 }

 //project initial precursors
 int cr_size = d_precursors_basis->number_columns();
 d_P0_r = new callow::Vector (cr_size, 0.0);
 d_precursors_basis->multiply_transpose(*d_P0, *d_P0_r);
}


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

 solver->solve(d_P_0, *d_P_r);
}


void TransientSolver::Solve(SP_state initial_state)
{
 d_material->update(0.0, 0, 1, false);
 d_state = initial_state;

 ProjectInitial(d_state);

 double t = 0.0;
 double dt = 0.0;

 int d_number_Steps; // TODO

 for (int i; i<-d_number_Steps; i++)
 {
  t += d_dt;
  step(t, dt);
  // update precursors

  // cyle states_precursors
 }

}


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
 SP_gainoperator G (new DiffusionGainOperator(d_inp, d_material, d_mesh, false));

 d_L = new callow::MatrixDense(rf, rf);

 OperatorProjection Projector(1);
 Projector.SetOperators(L, d_flux_basis);
 Projector.Project(d_L);

 d_G = new callow::MatrixDense(rf, rf);
 Projector.SetOperators(G, d_flux_basis);
 Projector.Project(d_G);

 //compute (I + dt*Pr)
  for (int i=0; i<d_rc; i++)
  {
   for (int j=0; j<d_rc ; j++)
   {
 	(*d_P)(i, j) = (*d_P)(i, j)*dt;
   }
   // add the identity matrix
   d_P->insert(i, i, 1.0, 1);
  }

  d_solver = LinearSolverCreator::Create(d_inp);

  d_solver->set_operators(d_P);
  // need to compute P^-1 * F -> PF
  MatrixDense PF(d_rc, d_rf);
  for (int i=0; i << d_rf; i++)
  {
   Vector y(d_rc, 0.0);
   Vector x(d_rc, 0.0);
   for (int j=0; j<< d_rc; j++)
   {
    x[j] = (*d_F)(j, i);
   }
   d_solver->solve(x, y);
   for (int k=0; k<d_rc; k++)
   {
    PF.insert(i, k, y[k]);
   }
  }

 // compute D@PF -> B
 d_B = new callow::Matrix (rf, rf);
 for (int i=0; i<rf ;i++)
 {
  for (int j=0; j<rc; j++)
  {
   double value = (*d_D)(i, j)*PF(j,i)*dt  ;
   d_B->insert(i, j, value, 1);
  }
 }

 d_A = new callow::MatrixDense(rf, rf);
 // compute the matrix
 for (int i=0; i<d_rf; i++)
 {
  for (int j=0; j<d_rf; j++)
  {
   double value = (*d_L)(i, j) - (*d_G)(i, j) - d_material->beta_total(0)*(*d_F)(i, j);
   d_A->insert(i, j, value);
  }
 }
}

void TransientSolver::compute_b()
{
 // compute P^-1*c0_r
 Vector Pc0(d_rc, 0.0);
 // need to change d_P0
 d_solver->solve(Pc0, *d_P0_r);

 // compute DP^-1*c0
 Vector DPc(d_rf, 0.0);
//DPc =  new Vector(rf, 0.0);
 d_D->multiply(Pc0, DPc);

 // compute b
 Vector b(d_rf, 0.0);

 for (int i=0; i<d_rf; i++)
 {
  b[i] = DPc[i] + (*d_phi0_r)[i];
 }

}

void TransientSolver::step(double t, double dt)
{
 // need to update the material here
 d_material->update(t, dt, 1, false, true);
 UpdateOperator(dt);

 d_solver->set_operators(d_A);
 d_solver->solve(*d_b, *d_phi_r);
 UpdatePrecursors(dt);

 // project back and store in d_phi if the final step
 if (t == 0)
 {
  d_flux_basis->multiply(*d_phi_r, *d_phi);
  d_precursors_basis->multiply(*d_P_r, *d_P0);
 }

 // cycle
 *d_phi0_r = *d_phi_r;
 *d_P0_r = *d_P_r;

}


