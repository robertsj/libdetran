//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGDiffusionSolver.cc
 *  @brief MGDiffusionSolver member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGDiffusionSolver.hh"
#include "boundary/BoundaryTraits.hh"
#include "DiffusionGainOperator.hh"
#include "callow/preconditioner/PCILU0.hh"
#include <cstdio>
#include <cmath>
#include <iostream>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
MGDiffusionSolver<D>::MGDiffusionSolver(SP_state                  state,
                                        SP_material               material,
                                        SP_boundary               boundary,
                                        const vec_externalsource &q_e,
                                        SP_fissionsource          q_f,
                                        bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
  , d_solver_type("gmres")
  , d_keff(1.0)
  , d_fill_boundary(false)
{
  // Set the problem dimension
  d_problem_size = d_mesh->number_cells() * d_material->number_groups();

  // Create vectors
  d_phi = new Vector_T(d_problem_size, 0.0);
  d_Q   = new Vector_T(d_problem_size, 0.0);

  // Create multigroup diffusion operator.  Note, the full energy range
  // is included.
  size_t cutoff = d_adjoint ? d_number_groups - 1 : 0;
  d_M   = new DiffusionLossOperator(d_input,
                                    d_material,
                                    d_mesh,
                                    d_multiply,
                                    cutoff,
                                    d_adjoint,
                                    1.0);       // Default to k=1

  // Get or create solver database.
  SP_input db;
  if (d_input->check("outer_solver_db"))
  {
    db = d_input->template get<SP_input>("outer_solver_db");
  }
  else
  {
    db = new detran_utilities::InputDB("mgdiffusionsolver_db");
    db->template put<double>("linear_solver_rtol", d_tolerance);
    db->template put<double>("linear_solver_atol", d_tolerance);
    db->template put<int>("linear_solver_maxit", d_maximum_iterations);
    db->template put<int>("linear_solver_monitor_level", d_print_level);
    d_input->template put<SP_input>("outer_solver_db", db);
  }

  // Build solver
  d_solver = Creator_T::Create(db);
  d_solver->set_operators(d_M, db);

  // Check whether we need boundary currents
  if (d_input->check("compute_boundary_flux"))
  {
    d_fill_boundary = d_input->template get<int>("compute_boundary_flux");
  }
}

//----------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::refresh()
{
  // Easiest approach for now is simply to rebuild the operator.  This
  // may be less than ideal, but the cost of building a diffusion operator
  // should almost always be a small fraction of the solver cost.
  d_M = new DiffusionLossOperator(d_input, d_material, d_mesh,
                                  d_multiply, 0, d_adjoint, d_keff);
  d_solver->set_operators(d_M);
}

//----------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::solve(const double keff)
{
  // Note, if keff is new, we must refresh the operators
  if (d_keff != keff)
  {
    d_keff = keff;
    refresh();
  }

  // Reset and build the right hand side.
  d_Q->set(0.0);
  build_volume_source();
  build_boundary_source();

  // Solve the problem
//  d_M->print_matlab("tran_diff.out");
//  d_M->compute_explicit("tran_exp.out");
  d_solver->solve(*d_Q, *d_phi);
//  d_Q->print_matlab("Q.out");
//  d_phi->print_matlab("phi.out");
//
//  d_Q->print_matlab("Q.out");
//  d_phi->print_matlab("phi.out");
//  d_M->print_matlab("M.out");

  // Fill the state and boundary
  fill_state();
  if (d_fill_boundary) fill_boundary();

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::build_volume_source()
{

  // Loop over external sources
  for (int i = 0; i < d_externalsources.size(); ++i)
  {
    if (!d_externalsources[i]) continue;
    // Loop over all groups.
    for (int g = 0; g < d_material->number_groups(); g++)
    {
      // Loop over all cells.
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        // Compute row index.
        int row = cell + g * d_mesh->number_cells();

        // Add source
        (*d_Q)[row] += d_externalsources[i]->source(cell, g);

      } // row loop
    } // group loop
  } // source loop

  // Add fission source, if applicable. This only happens for
  // eigenvalue problems or multiplying fixed source problems
  // in which fission iterations are performed.  In both cases,
  // the caller has updated the source.
  if (d_fissionsource && !d_multiply)
  {
    for (int g = 0; g < d_material->number_groups(); g++)
    {
      // Loop over all cells.
      for (int cell = 0; cell < d_mesh->number_cells(); cell++)
      {
        // Compute row index.
        int row = cell + g * d_mesh->number_cells();

        // Add fission source
        (*d_Q)[row] += d_fissionsource->source(g)[cell];

      } // row loop
    } // group loop
  }

}

//---------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::build_boundary_source()
{

  // Make reference to boundary for better notation
  typename Boundary_T::SP_boundary b_sp = d_boundary;
  Boundary_T &J = *b_sp;

  // For a given dimension, provide remaining dimensions
  int remdims[3][2] = {{1,2}, {0,2}, {0,1}};

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

//---------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::fill_state()
{
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
  if (d_state->store_current())
  {
    fill_current();
  }
}

//---------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::fill_current()
{
  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    moments_type &J_g   = d_state->current(g);
    moments_type &phi_g = d_state->phi(g);

    for (int cell = 0; cell < d_mesh->number_cells(); ++cell)
    {

      // Define the data for this cell.
      size_t m = mat_map[cell];

      double cell_dc = d_material->diff_coef(m, g);
      Assert(cell_dc > 0.0);

      // Get the directional indices.
      size_t i = d_mesh->cell_to_i(cell);
      size_t j = d_mesh->cell_to_j(cell);
      size_t k = d_mesh->cell_to_k(cell);

      // Cell width vector.
      double cell_hxyz[3] = {d_mesh->dx(i), d_mesh->dy(j), d_mesh->dz(k)};

      // Face currents.
      double J[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      // Index arrays to help determine if a cell surface is on the boundary.
      int bound[6] = {i, i, j, j, k, k};
      int nxyz[3][2] = {{0, d_mesh->number_cells_x()-1},
                        {0, d_mesh->number_cells_y()-1},
                        {0, d_mesh->number_cells_z()-1}};

      // leak --> 0=-x, 1=+x, 2=-y, 3=+y, 4=-z, 5=+z
      for (int leak = 0; leak < 6; leak++)
      {
        int xyz_idx = std::floor(leak / 2); // determine W/E, S/N, B/T
        int dir_idx = 1 - ((leak + 1) % 2); // e.g. W vs E

        // Determine whether the neighbor is positive (+1) or negative (-1)
        // relative to the surface under consideration, and then decrement
        // the appropriate x, y, or z index.
        int neig_idx[3] = {i, j, k};
        int shift_idx   = -2 * ((leak + 1) % 2) + 1;
        neig_idx[xyz_idx] += shift_idx;

        // Compute coupling coefficient
        double dtilde = 0.0;
        if (bound[leak] == nxyz[xyz_idx][dir_idx])
        {

          dtilde = ( 2.0 * cell_dc * (1.0 - d_M->albedo(leak, g)  ) ) /
                   ( 4.0 * cell_dc * (1.0 + d_M->albedo(leak, g)) +
                    (1.0 - d_M->albedo(leak, g)) * cell_hxyz[xyz_idx] );

          J[leak] = dtilde * phi_g[cell];

        }
        else // not a boundary
        {
          // Get the neighbor data.
          size_t neig_cell = d_mesh->index(neig_idx[0], neig_idx[1], neig_idx[2]);
          size_t ii = d_mesh->cell_to_i(neig_cell);
          size_t jj = d_mesh->cell_to_j(neig_cell);
          size_t kk = d_mesh->cell_to_k(neig_cell);

          // Neighbor volume and diffusion coefficient.
          double neig_hxyz[3] = {d_mesh->dx(ii), d_mesh->dy(jj), d_mesh->dz(kk)};
          double neig_dc = d_material->diff_coef(mat_map[neig_cell], g);

          // Compute dtilde.
          dtilde = ( 2.0 * cell_dc * neig_dc ) /
                   ( neig_hxyz[xyz_idx] * cell_dc +
                     cell_hxyz[xyz_idx] * neig_dc );

          // Compute surface current
          J[leak] = (double)dir_idx * dtilde * (phi_g[neig_cell] - phi_g[cell]);

        }

      } // leak loop

      double A = cell_hxyz[1] * cell_hxyz[2] +
                 cell_hxyz[0] * cell_hxyz[2] +
                 cell_hxyz[0] * cell_hxyz[1];
      A *= 1.0/6.0;
      double J_x = (J[0] + J[1]) * cell_hxyz[1] * cell_hxyz[2] / A;
      double J_y = (J[2] + J[3]) * cell_hxyz[0] * cell_hxyz[2] / A;
      double J_z = (J[4] + J[5]) * cell_hxyz[0] * cell_hxyz[1] / A;

      // Area-averaged current norm
      J_g[cell] = std::sqrt(J_x*J_x + J_y*J_y + J_z*J_z);

    }

  }

}

//---------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::fill_boundary()
{
  // Reference to boundary for better notation
  typename Boundary_T::SP_boundary b_sp = d_boundary;
  Boundary_T &J = *b_sp;

  // For a given dimension, provide remaining dimensions
  int remdims[3][2] = {{1,2}, {0,2}, {0,1}};

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
              ((2.0 * DC) * (*d_phi)[row] + (W - 4.0*DC) * Jinc) /
                (4.0 * DC + W);

          } // end group loop
        } // end dim2 loop
      } // end dim1 loop
    } // end dir loop
  } // end dim0 loop
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class MGDiffusionSolver<_1D>;
template class MGDiffusionSolver<_2D>;
template class MGDiffusionSolver<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGDiffusionSolver.cc
//----------------------------------------------------------------------------//
