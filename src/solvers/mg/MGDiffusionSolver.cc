//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGDiffusionSolver.cc
 *  @brief  MGDiffusionSolver member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 11, 2012
 */
//---------------------------------------------------------------------------//

#include "MGDiffusionSolver.hh"
#include "boundary/BoundaryTraits.hh"
#include "DiffusionGainOperator.hh"
#include "callow/preconditioner/PCILU0.hh"
#include <cstdio>
#include <cmath>
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
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
{

  // Set the problem dimension
  d_problem_size = d_mesh->number_cells() * d_material->number_groups();

  // Create vectors
  d_phi = new Vector_T(d_problem_size, 0.0);
  d_Q   = new Vector_T(d_problem_size, 0.0);

  // Create multigroup diffusion operator
  d_M   = new DiffusionLossOperator(d_input,
                                    d_material,
                                    d_mesh,
                                    d_multiply,
                                    0,          // Full energy spectrum
                                    d_adjoint,
                                    1.0);       // Default to k=1

  // Create solver and set operator.
  SP_input db;
  if (d_input->check("outer_solver_db"))
  {
    db = d_input->template get<SP_input>("outer_solver_db");
  }
  d_solver = Creator_T::Create(db);
  d_solver->set_operators(d_M, db);

}

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
template <class D>
void MGDiffusionSolver<D>::solve(const double keff)
{
  // Note, if keff is new, the client is responsible to refresh
  // the solver.
  d_keff = keff;

  // Reset and build the right hand side.
  d_Q->set(0.0);
  build_volume_source();
  build_boundary_source();

  // Solve the problem
  d_solver->solve(*d_Q, *d_phi);
//  d_M->print_matlab("M.out");
//  d_phi->print_matlab("phi.out");
//  d_M->display();
//  d_Q->display();
//  d_phi->display();
//  d_M->multiply(*d_phi, *d_Q);
//  d_Q->display();

//  THROW("done");
  // Fill the state and boundary
  fill_state();
  fill_boundary();

  //d_state->display();
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
  if (d_fissionsource and !d_multiply)
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

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class MGDiffusionSolver<_1D>;
template class MGDiffusionSolver<_2D>;
template class MGDiffusionSolver<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file MGDiffusionSolver.cc
//---------------------------------------------------------------------------//
