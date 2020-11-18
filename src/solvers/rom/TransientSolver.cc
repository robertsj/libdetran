//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TransientSolver.cc
 *  @brief TransientSolver class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


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

  d_P0 = new callow::Vector (d_num_cells* d_precursors_group, 0.0);
  d_phi0 = new callow::Vector(d_num_cells* d_number_groups, 0.0);

  d_phi0_r = new callow::Vector(d_rf, 0.0);
  d_P0_r = new callow::Vector (d_rc, 0.0);

  // matrix of precursors concentraion at all time steps
  d_precursors = new callow::MatrixDense(d_num_cells*d_precursors_group, d_number_steps+1);
  // matrix of the flux at all time steps
  d_flux = new callow::MatrixDense(d_number_groups*d_num_cells, d_number_steps+1);

  d_A = new callow::MatrixDense(d_rf + d_rc, d_rf+d_rc);

  // long vector of flux and precursors
  d_sols = new callow::MatrixDense(d_rf + d_rc, d_number_steps+1);

  // long vector of reduced flux and precursors
  d_sol_r = new callow::Vector(d_rf+d_rc, 0.0);

  d_sol0_r = new callow::Vector(d_rf + d_rc, 0.0);

  if (d_inp->check("rom_solver_db"))
  {
    db = d_inp->template get<SP_input>("rom_solver_db");
  }

  d_solver = LinearSolverCreator::Create(db);
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
     // store the initial concentration in the solution matrix
     (*d_precursors)(cell + i*d_num_cells, 0) =  (*d_P0)[cell + i*d_num_cells];
   }
 }
}

//------------------------------------------------------------------------------------//

/// Project the initial flux and precursors concentraion onto their basis to get
void TransientSolver::ProjectInitial()
{
  for (int g=0; g<d_number_groups; g++)
  {
    for (int cell=0; cell<d_num_cells; cell++)
    {
      (*d_phi0)[cell + g*d_num_cells] = d_state->phi(g)[cell];
      //store the initial flux
      (*d_flux)(cell + g*d_num_cells, 0) = d_state->phi(g)[cell];
    }
  }

  d_flux_basis->multiply_transpose(*d_phi0, *d_phi0_r);

  // project initial precursors
  d_precursors_basis->multiply_transpose(*d_P0, *d_P0_r);

  // stack the flux and precursors in one vector
  for (int i =0; i<d_rf; i++)
  {
    (*d_sol0_r)[i] = (*d_phi0_r)[i];
  }

  for (int i =0; i<d_rc; i++)
  {
    (*d_sol0_r)[d_rf + i] = (*d_P0_r)[i];
  }
}

//------------------------------------------------------------------------------------//

void TransientSolver::Construct_Operator(double t, double dt)
{
  // get the matrices
  KineticMatrices K(d_inp, d_mesh, d_material, d_flux_basis, d_precursors_basis);

  d_precursors_decay = new callow::MatrixDense(d_rc, d_rc);
  d_precursors_decay = K.precursors_decay();

  d_delayed_production = new callow::MatrixDense(d_rf, d_rc);
  d_delayed_production = K.delayed_production();

  d_precursors_production = new callow::MatrixDense(d_rc, d_rf);
  d_precursors_production = K.precursors_production();

  // if diffusion
  SP_lossoperator L(new DiffusionLossOperator(d_inp, d_material, d_mesh, false, 0.0, false, 1.0));
  d_L = L;

  SP_gainoperator G(new DiffusionGainOperator(d_inp, d_material, d_mesh, false));
  d_G = G;
  //d_G = new DiffusionGainOperator(d_inp, d_material, d_mesh, false);

 callow::Matrix::SP_matrix d_L_prime;
 callow::Matrix::SP_matrix d_G_prime;

  int m =d_number_groups*d_num_cells;

  vec_int nnz_L(m, 1 + 2 * d_mesh->dimension() + d_number_groups);
  d_L_prime = new callow::Matrix(d_number_groups*d_num_cells, d_number_groups*d_num_cells, nnz_L[0]);

  vec_int nnz_G(m, d_number_groups);
  d_G_prime = new callow::Matrix(d_number_groups*d_num_cells, d_number_groups*d_num_cells, nnz_G[0]);

  int * rows_L = d_L->rows();
  int* cols_L = d_L->columns();
  double* v_L = d_L->values();

  int * rows_G = d_G->rows();
  int* cols_G = d_G->columns();
  double* v_G = d_G->values();

  OperatorProjection Projector(1);
  d_Gr = new callow::MatrixDense(d_rf, d_rf);
  d_Lr = new callow::MatrixDense(d_rf, d_rf);


  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

  for (int g=0; g < d_number_groups; g++)
  {
    // loop over cells
    for (int cell=0; cell < d_num_cells; cell++)
    {
      int m = mat_map[cell];
      int r = d_num_cells*g + cell; //row number
      for (int p = rows_L[r]; p < rows_L[r + 1]; p++)
      {
	int c = cols_L[p];
        double value = -v_L[p]*d_material->velocity(g);
        //d_L->insert(r, c, value, 0);
        d_L_prime->insert(r, c, value);
      }

      for (int p = rows_G[r]; p < rows_G[r + 1]; p++)
      {
        int c = cols_G[p];
        double value = v_G[p]*(1 - d_material->beta_total(m))*d_material->velocity(g);
        //d_G->insert(r, c, value, 0);
        d_G_prime->insert(r, c, value);
      }
    }
 }

 d_L_prime->assemble();
 d_G_prime->assemble();

  Projector.SetOperators(d_G_prime, d_flux_basis);

  Projector.Project(d_Gr);


  Projector.SetOperators(d_L_prime, d_flux_basis);
  Projector.Project(d_Lr);



  // this assumes that a basis set is generated for each flux group, so the velocity
  // is not collapsed.
  int r = d_rf/d_number_groups;
  for (int g=0; g<d_number_groups; g++)
  {
    for (int i=0; i<r; i++)
    { // fill the upper left of A
      for (int j=0; j<d_rf; j++)
      {
        (*d_A)(i + g*r, j) = ((*d_Lr)(i + g*r, j) + (*d_Gr)(i + g*r, j))*-d_dt;
         if (i + g*r == j) (*d_A)(i + g*r, j) += 1;
      }
      // fill the upper right
      for (int k =0; k<d_rc; k++)
      {
        (*d_A)(i + g*r, k+d_rf) = (*d_delayed_production)(i+g*r, k)*d_material->velocity(g)*-d_dt;
        if (i + g*r == k+d_rf) (*d_A)(i + g*r, k+d_rf) += 1;
      }
    }
  }
   // fill the lower half
  for (int i=0; i<d_rc ; i++)
  {  // lower left
     for (int j=0; j<d_rf; j++)
     {
        (*d_A)(d_rf+i, j) = (*d_precursors_production)(i, j)*-d_dt;
        if (d_rf+i == j) (*d_A)(d_rf+i, j) += 1;
     }
    // lower right
    for (int k=0; k<d_rc; k++)
    {
      (*d_A)(d_rf+i, k+d_rf) = (*d_precursors_decay)(i, k)*-d_dt;
      if (d_rf+i == k+d_rf) (*d_A)(d_rf+i, k+d_rf) += 1;
    }
  }
 }

//---------------------------------------------------------------------//

void TransientSolver::Refersh_Operator()
{
  // if diffusion
  SP_lossoperator L(new DiffusionLossOperator(d_inp, d_material, d_mesh, false, 0.0, false, 1.0));
  d_L = L;

  int * rows_L = d_L->rows();
  int* cols_L = d_L->columns();
  double* v_L = d_L->values();

  d_Lr = new callow::MatrixDense(d_rf, d_rf);

  OperatorProjection Projector(1);
  Projector.SetOperators(d_L, d_flux_basis);
  Projector.Project(d_Lr);
  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

  // this assumes that a basis set is generated for each flux group, so the velocity
  // is not collapsed.
  int r = d_rf/d_number_groups;

  for (int g=0; g<d_number_groups; g++)
  {
    for (int i=0; i<r; i++)
    {
      for (int j=0; j<d_rf; j++)
      {
        (*d_A)(i + g*r, j) = (-(*d_Lr)(i + g*r, j)*d_material->velocity(g) + (*d_Gr)(i + g*r, j))*-d_dt;
        if (i + g*r == j) (*d_A)(i + g*r, j) += 1;
      }
    }
  }
}
//------------------------------------------------------------------------------------//

void TransientSolver::Solve(SP_state initial_state)
{
  d_material->update(0.0, 0, 1, false);

  Construct_Operator(0, d_dt);

  d_state = initial_state;

  initialize_precursors();

  ProjectInitial();

  double t = 0.0;
  int n = d_rf + d_rc;

  SP_matrix d_A_;
  d_A_ = new callow::MatrixDense(n, n);

  for (int step=0 ; step< d_number_steps; step++)
  {
    t += d_dt;
    d_material->update(t, d_dt, 1, false);

    if (d_material->update_material)
    {
     Refersh_Operator();
    }

   d_solver->set_operators(d_A, db);
   d_solver->solve(*d_sol0_r, *d_sol_r);

   // store this state in the solution matrix
   for (int i=0; i<n; i++)
   {
     (*d_sols)(i, step) = (*d_sol_r)[i];
     (*d_sol0_r)[i] = (*d_sol_r)[i];
   }
 }

 // reconstruct the full solution
  reconstruct();
}

//------------------------------------------------------------------------------------//

void TransientSolver::reconstruct()
{
  callow::Vector phi(d_number_groups*d_num_cells, 0.0);
  callow::Vector C(d_precursors_group*d_num_cells, 0.0);

  callow::Vector v1(d_rf, 0.0);
  callow::Vector v2(d_rc, 0.0);

  for (int i=0; i<d_number_steps; i++)
  {
    for (int f=0; f< d_rf; f++)
    {
      v1[f] = (*d_sols)(f, i);
    }
    d_flux_basis->multiply(v1, phi);

    for (int p=0; p<d_rc; p++)
    {
      v2[p] = (*d_sols)(p + d_rf, i);
    }

    d_precursors_basis->multiply(v2, C);
    double *phi_ = &phi[0];
    double *C_ = &C[0];

    d_flux->insert_col(i+1, phi_, 0);
    d_precursors->insert_col(i+1, C_, 0);
  }
  // temporary ... need to have getter for flux, etc
  d_flux->print_matlab("flux.txt");
}
