//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionFixedSourceSolve.cc
 * \brief  DiffusionFixedSourceSolve member definitions
 * \author Jeremy Roberts
 * \date   Sep 11, 2012
 */
//---------------------------------------------------------------------------//

#include "DiffusionFixedSourceSolver.hh"
#include "boundary/BoundaryTraits.hh"
#include <cstdio>
#include <cmath>
#include <iostream>

namespace detran
{

template <class D>
DiffusionFixedSourceSolver<D>::DiffusionFixedSourceSolver(SP_input input,
                                                          SP_material material,
                                                          SP_mesh mesh,
                                                          SP_state state)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_state(state)
  , d_boundary(new BoundaryDiffusion<D>(input, mesh))
  , d_fixed_type(FIXED)
  , d_maximum_iterations(100)
  , d_maximum_fission_iterations(20)
  , d_tolerance(1e-6)
  , d_fission_tolerance(1e-6)
  , d_fission_scaling(1.0)
  , d_print_out(0)
  , d_print_interval(10)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_state);
  Require(d_boundary);

  d_problem_size = d_mesh->number_cells() * d_material->number_groups();

  // Select solver type
  if (d_input->check("diffusion_fixed_type"))
  {
    d_fixed_type = d_input->template get<int>("diffusion_fixed_type");
    Ensure(d_fixed_type < END_FIXED_TYPES);
  }

  // Create operators and vectors
  d_phi = new Vector(d_problem_size, 0.0);
  d_Q   = new Vector(d_problem_size, 0.0);
  if (d_fixed_type == MULTIPLY)
    d_M = new DiffusionLossOperator(d_input, d_material, d_mesh, true, d_fission_scaling);
  else
    d_M =  new DiffusionLossOperator(d_input, d_material, d_mesh, false);
  if (d_fixed_type == ITERATE)
  {
    d_F = new DiffusionGainOperator(d_input, d_material, d_mesh);
    d_phi_old = new Vector(d_problem_size, 0.0);
    d_Q_total = new Vector(d_problem_size, 0.0);
  }

  // Create solver
  PetscErrorCode ierr;
  ierr = KSPCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating KSP object.");

  // Set the operators.  We use the same operator for the PC.
  KSPSetOperators(d_solver, d_M->A(), d_M->A(), SAME_NONZERO_PATTERN);

  // Set tolerances.
  ierr = KSPSetTolerances(d_solver,
                          d_tolerance,   // relative tolerance
                          PETSC_DEFAULT, // absolute tolerance
                          PETSC_DEFAULT, // divergence tolerance
                          d_maximum_iterations);

  // Allow for command line flags.
  ierr = KSPSetFromOptions(d_solver);

  d_M->display();

  // Postconditions
  Ensure(!ierr);
}

template <class D>
void DiffusionFixedSourceSolver<D>::solve()
{
  // Call the appropriate solver
  if (d_fixed_type == FIXED or d_fixed_type == MULTIPLY)
    solve_fixed();
  else if (d_fixed_type == ITERATE)
    solve_iterate();
  else
    THROW("Unsupported diffusion fixed source type");

  std::cout << " filling state..." << std::endl;
  // Fill the state with the converged flux
  size_t row = 0;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    moments_type &phi_g = d_state->phi(g);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell, ++row)
    {
      phi_g[cell] = (*d_phi)[row];
    }
  }

  // Fill the boundary with the outgoing current
  std::cout << " filling boundary..." << std::endl;
  fill_boundary();
}

template <class D>
void DiffusionFixedSourceSolver<D>::build_source(SP_source q)
{
  if (q)
  {
    build_volume_source(q);
    std::cout << "built volume source" << std::endl;
  }
  build_boundary_source();
  std::cout << "built boundary source" << std::endl;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

// Direct solve
template <class D>
void DiffusionFixedSourceSolver<D>::solve_fixed()
{
  PetscErrorCode ierr = KSPSolve(d_solver, d_Q->V(), d_phi->V());
  Insist(!ierr, "Error in KSPSolve.");
}

// Fission source iterations
template <class D>
void DiffusionFixedSourceSolver<D>::solve_iterate()
{
  // Clear the total source and add the fixed source
  d_Q_total->scale(0.0);
  d_Q_total->add(*d_Q);

  // Pointer for swapping
  SP_vector phi_swap;

  // Norm of the residual
  double error = 0.0;

  for (int iteration = 0; iteration < d_maximum_fission_iterations; ++iteration)
  {

    //-----------------------------------------------------------------------//
    // SWAP OLD AND NEW
    //-----------------------------------------------------------------------//

    phi_swap  = d_phi;
    d_phi     = d_phi_old;
    d_phi_old = phi_swap;

    //-----------------------------------------------------------------------//
    // SOLVE phi(n+1) = inv(M) * b_total(n+1)
    //-----------------------------------------------------------------------//

    // Compute
    PetscErrorCode ierr = KSPSolve(d_solver, d_Q_total->V(), d_phi->V());
    Insist(!ierr, "Error in KSPSolve.");

    //-----------------------------------------------------------------------//
    // CHECK CONVERGENCE
    //-----------------------------------------------------------------------//

    error = d_phi->residual_norm(*d_phi_old, Vector::L2);
    if (d_print_out > 1 and iteration % d_print_interval == 0)
    {
      printf("  Fixed fission iteration: %3i  Error: %12.9f \n",
             iteration, error);
    }
    if (error < d_fission_tolerance) break;

    //-----------------------------------------------------------------------//
    // UPDATE FISSION SOURCE
    //-----------------------------------------------------------------------//

    // Q_total <-- (F*x) / k
    d_F->multiply(*d_phi, *d_Q_total);
    d_Q_total->scale(1/d_fission_scaling);
    // Q_total <-- (F*x) / k + Q
    d_Q_total->add(*d_Q);

  }

}

template <class D>
void DiffusionFixedSourceSolver<D>::build_volume_source(SP_source q)
{
  // Preconditions
  Require(q);
  std::cout << "building v source..." << std::endl;

  // Loop over all groups.
  for (int g = 0; g < d_material->number_groups(); g++)
  {
    // Loop over all cells.
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {

      // Compute row index.
      int row = cell + g * d_mesh->number_cells();

      // Add source
      (*d_Q)[row] += q->source(cell, g);

    } // row loop

  } // group loop
}

template <class D>
void DiffusionFixedSourceSolver<D>::build_boundary_source()
{
  std::cout << "building b source..." << std::endl;
  // Make reference to boundary for better notation
  Boundary_T &J = *d_boundary;

  // For a given dimension, provide remaining dimensions
  int remdims[3][2] = {{2,1}, {2,0}, {1,0}};

  // Cell indices
  int ijk[3] = {0, 0, 0};
  int &i = ijk[0];
  int &j = ijk[1];
  int &k = ijk[2];

  // Get the material map.
  detran_utilities::vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  // Loop over all dimensions
  for (int dim0 = 0; dim0 < D::dimension; ++dim0)
  {
    // Bounding cell indices for this dimension
    int bound[2] = {0, d_mesh->number_cells(dim0)-1};

    // Other dimensions
    int dim1 = remdims[dim0][0];
    int dim2 = remdims[dim0][1];

    // Loop over directions - and +
    for (int dir = 0; dir < 2; ++dir)
    {
      // Surface index
      int surface = 2 * dim0 + dir;

      // Index and width along this direction
      ijk[dim0] = bound[dir];
      double W  = d_mesh->width(dim0, ijk[dim0]);

      // Loop over secondaries dimensions on this surface
      for (ijk[dim1] = 0; ijk[dim1] < d_mesh->number_cells(dim1); ++ijk[dim1])
      {
        for (ijk[dim2] = 0; ijk[dim2] < d_mesh->number_cells(dim2); ++ijk[dim2])
        {
          // Get the cardinal cell index
          size_t cell = d_mesh->index(i, j, k);

          // Material index for this cell
          size_t m = mat_map[cell];

          // Loop over all groups
          for (int g = 0; g < d_material->number_groups(); g++)
          {
            // Compute row index of system
            int row = cell + g * d_mesh->number_cells();

            // Diffusion coefficient
            double DC = d_material->diff_coef(m, g);

            // Incident partial current
            double Jinc = BoundaryValue<D>::value
                          (J(surface, g, Boundary_T::IN), ijk[dim1], ijk[dim2]);

            // Compute source for this cell
            (*d_Q)[row] += (8.0 * DC * Jinc) / ((4.0 * DC + W) * W);

          } // end group loop

        } // end dim2 loop

      } // end dim1 loop

    } // end dir loop

  } // end dim0 loop

}


template <class D>
void DiffusionFixedSourceSolver<D>::fill_boundary()
{
  // Reference to boundary for better notation
  Boundary_T &J = *d_boundary;

  // For a given dimension, provide remaining dimensions
  int remdims[3][2] = {{2,1}, {2,0}, {1,0}};

  // Cell indices
  int ijk[3] = {0, 0, 0};
  int &i = ijk[0];
  int &j = ijk[1];
  int &k = ijk[2];

  // Get the material map.
  detran_utilities::vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  // Loop over all dimensions
  for (int dim0 = 0; dim0 < D::dimension; ++dim0)
  {
    // Bounding cell indices for this dimension
    int bound[2] = {0, d_mesh->number_cells(dim0)-1};

    // Other dimensions
    int dim1 = remdims[dim0][0];
    int dim2 = remdims[dim0][1];

    // Loop over directions - and +
    for (int dir = 0; dir < 2; ++dir)
    {
      // Surface index
      int surface = 2 * dim0 + dir;

      // Index and width along this direction
      ijk[dim0] = bound[dir];
      double W  = d_mesh->width(dim0, ijk[dim0]);

      // Loop over secondaries dimensions on this surface
      for (ijk[dim1] = 0; ijk[dim1] < d_mesh->number_cells(dim1); ++ijk[dim1])
      {
        for (ijk[dim2] = 0; ijk[dim2] < d_mesh->number_cells(dim2); ++ijk[dim2])
        {
          // Get the cardinal cell index
          size_t cell = d_mesh->index(i, j, k);

          // Material index for this cell
          size_t m = mat_map[cell];

          // Loop over all groups
          for (int g = 0; g < d_material->number_groups(); g++)
          {
            // Compute row index of system
            int row = cell + g * d_mesh->number_cells();

            // Diffusion coefficient
            double DC = d_material->diff_coef(m, g);

            // Incident partial current
            double Jinc = BoundaryValue<D>::value
                          (J(surface, g, Boundary_T::IN), ijk[dim1], ijk[dim2]);

            // Compute the outgoing partial current
            BoundaryValue<D>::value
            (J(surface, g, Boundary_T::OUT), ijk[dim1], ijk[dim2]) =
              ((2.0 * DC) * (*d_phi)[row] + (W - 4.0*DC) * Jinc) / (4.0 * DC + W);

          } // end group loop

        } // end dim2 loop

      } // end dim1 loop

    } // end dir loop

  } // end dim0 loop
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class DiffusionFixedSourceSolver<_1D>;
template class DiffusionFixedSourceSolver<_2D>;
template class DiffusionFixedSourceSolver<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file DiffusionFixedSourceSolve.cc
//---------------------------------------------------------------------------//
