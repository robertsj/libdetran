//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CMFDLossOperator.cc
 *  @brief CMFDLossOperator
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "CMFDLossOperator.hh"
#include "utilities/MathUtilities.hh"
#include <iostream>

#define COUT(c) std::cout << c << std::endl;

namespace detran
{

using std::cout;
using std::endl;
using detran_utilities::range;

//----------------------------------------------------------------------------//
template <class D>
CMFDLossOperator<D>::CMFDLossOperator(SP_input      input,
                                      SP_material   material,
                                      SP_mesh       mesh,
                                      SP_tally      tally,
                                      const bool    include_fission,
                                      const bool    adjoint,
                                      const double  keff)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_tally(tally)
  , d_include_fission(include_fission)
  , d_adjoint(adjoint)
  , d_keff(keff)
  , d_correct(true)
  , d_initial(true)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_tally);

  // Set the dimensions
  d_dimension = d_mesh->dimension();
  d_number_groups = d_material->number_groups();
  int upper = d_adjoint ? -1 : d_number_groups;
  int lower = d_adjoint ? d_number_groups - 1 : 0;
  d_groups = detran_utilities::range<size_t>(lower, upper);
  d_group_size = d_mesh->number_cells();
  Base::set_size(d_number_groups * d_group_size,
                 d_number_groups * d_group_size);

  // Default albedos to 1.0.  For dimensions in play, this will be
  // overwritten by the default boundary.
  d_albedo.resize(6,  vec_dbl(d_number_groups, 1.0));

  // Preallocate.  The number of nonzeros is
  //   diagonal + 2*dim neighbors + num_groups coupling from scatter/fission
  vec_int nnz(d_m, 1 + 2 * d_dimension + d_number_groups);
  preallocate(&nnz[0]);

  // Set the albedo.  First, check if the input has an albedo
  // entry.  If it does, this is the default way to set the
  // condition.  Otherwise, check the boundary conditions.
  bool zero_flux = false;
  if (d_input->check("bc_zero_flux"))
    zero_flux = d_input->get<int>("bc_zero_flux");
  if (d_input->check("albedo"))
  {
    // \todo Add a user-defined albedo
  }
  else
  {
    std::vector<std::string> boundary_name(6, "");
    boundary_name[Mesh::WEST]   = "bc_west";
    boundary_name[Mesh::EAST]   = "bc_east";
    boundary_name[Mesh::SOUTH]  = "bc_south";
    boundary_name[Mesh::NORTH]  = "bc_north";
    boundary_name[Mesh::BOTTOM] = "bc_bottom";
    boundary_name[Mesh::TOP]    = "bc_top";
    for (groups_iter g = d_groups.begin(); g != d_groups.end(); ++g)
    {
      for (int b = 0; b < d_mesh->dimension() * 2; b++)
      {
        d_albedo[b][*g] = 0.0;
        if (d_input->check(boundary_name[b]))
        {
          if (d_input->get<std::string>(boundary_name[b]) == "reflect")
          {
            d_albedo[b][*g] = 1.0;
          }
          else
          {
            if (zero_flux) d_albedo[b][*g] = -1.0;
          }
        }
      }
    }
  }
  // note, needs to be constructed

  if (d_input->check("cmfd_correct_coupling"))
    d_correct = 0 != d_input->get<int>("cmfd_correct_coupling");

  if (d_input->check("cmfd_relaxation"))
    d_alpha = d_input->get<double>("cmfd_relaxation");
  Require(d_alpha > 0.0 && d_alpha <= 1.0);

  // Coupling coefficient
  d_d_hat.resize(d_number_groups,
                 vec2_dbl(d_group_size,
                          vec_dbl(d_dimension * 2, 0.0)));
}

//----------------------------------------------------------------------------//
template <class D>

void CMFDLossOperator<D>::construct(const vec2_dbl &phi,
                                    const double    keff,
                                    SP_material     mat,
                                    bool            init)
{
  Require(phi.size() == d_number_groups);
  Require(phi[0].size() == d_group_size);
  clear();
  d_keff = keff;
  if (mat)
  {
    Require(mat->number_groups() == d_material->number_groups());
    Require(mat->number_materials() == d_material->number_materials());
    d_material = mat;
  }
  d_initial = init;
  build(phi);
}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
template <class D>
void CMFDLossOperator<D>::build(const vec2_dbl &phi)
{
  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");

  double alpha = d_alpha;
  if (d_initial) alpha = 1.0;

  bool corrected = false;
  for (groups_iter g_it = d_groups.begin(); g_it != d_groups.end(); ++g_it)
  {
    int g = *g_it;

    // Loop over all cells.
    for (int cell = 0; cell < d_group_size; cell++)
    {
      // Compute row index.
      int row = cell + g * d_group_size;

      // Define the data for this cell.
      size_t m = mat_map[cell];
      double cell_dc = d_material->diff_coef(m, g);
      Assert(cell_dc != 0.0);
      double cell_sr = d_material->sigma_t(m, g);
      cell_sr -= d_material->sigma_s(m, g, g);

      // Get the directional indices.
      size_t i = d_mesh->cell_to_i(cell);
      size_t j = d_mesh->cell_to_j(cell);
      size_t k = d_mesh->cell_to_k(cell);

      // Cell width vector.
      double cell_hxyz[3] = {d_mesh->dx(i), d_mesh->dy(j), d_mesh->dz(k)};

      // Direction-specific leakage coefficients.
      double jo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      // Index arrays to help determine if a cell surface is on the boundary.
      int bound[6] = {i, i, j, j, k, k};
      int nxyz[3][2] = {{0, d_mesh->number_cells_x()-1},
                        {0, d_mesh->number_cells_y()-1},
                        {0, d_mesh->number_cells_z()-1}};

      // Flag for insertion status.  If false, it means we didn't add
      // the value and should counter accordingly (here, we assert)
      bool flag;

      // For each spatial cell, there are 6 faces that connect the
      // cell to a neighbor or the global boundary.  Looping through
      // each of the faces, we can use the indices of the cell to
      // determine whether the face is between cells or is a boundary
      // surface.

      // leak --> 0=-x, 1=+x, 2=-y, 3=+y, 4=-z, 5=+z
      for (int leak = 0; leak < D::dimension*2; ++leak)
      {
        // Determine whether this is a W/E, S/N, or B/T face
        int xyz_idx = std::floor(leak / 2);

        // Determine the direction, e.g. left (-) vs right (+).
        int dir_idx = 1 - ((leak + 1) % 2);

        // Put i, j, k into an array.  The neighbor indices are
        // a perturbation of this.
        int neig_idx[3] = {i, j, k};
        int &ii = neig_idx[0];
        int &jj = neig_idx[1];
        int &kk = neig_idx[2];
        int ijk = neig_idx[xyz_idx];

        // Determine whether the neighbor is positive (+1) or negative (-1)
        // relative to the surface under consideration, and then decrement
        // the appropriate x, y, or z index.
        int shift_idx   = -2 * ((leak + 1) % 2) + 1;
        neig_idx[xyz_idx] += shift_idx;

        // Compute coupling coefficient
        double dtilde = 0.0;
        double dhat = 0.0;

        // net current leaving this cell through the indexed surface
        int surf_idx[] = {i, j, k};
        surf_idx[xyz_idx] += leak % 2;
        int si = surf_idx[0];
        int sj = surf_idx[1];
        int sk = surf_idx[2];
        double J = d_tally->partial_current(si, sj, sk, g, xyz_idx, true)
                 - d_tally->partial_current(si, sj, sk, g, xyz_idx, false);
        J *= (double)shift_idx * d_mesh->width(xyz_idx, ijk) / d_mesh->volume(cell);

        // cell flux
        double phi_cell = phi[g][cell];

        if (bound[leak] == nxyz[xyz_idx][dir_idx])
        {
          // Compute coupling coefficients
          dtilde = ( 2.0 * cell_dc * (1.0 - d_albedo[leak][g]) )  /
                   ( 4.0 * cell_dc * (1.0 + d_albedo[leak][g]) +
                    (1.0 - d_albedo[leak][g]) * cell_hxyz[xyz_idx]);
          dhat = (J - dtilde * phi_cell) / phi_cell;

          // Relax
          d_d_hat[g][cell][leak] = alpha * dhat
                                 + (1.0 - alpha) * d_d_hat[g][cell][leak];

        }
        else
        {
          // Get the neighbor data.
          size_t neig_cell = d_mesh->index(neig_idx[0], neig_idx[1], neig_idx[2]);

          // Neighbor volume, diffusion coefficient, and flux
          double neig_hxyz[3] = {d_mesh->dx(ii), d_mesh->dy(jj), d_mesh->dz(kk)};
          double neig_dc = d_material->diff_coef(mat_map[neig_cell], g);
          double phi_neig = phi[g][neig_cell];

          // Compute coupling coefficients
          dtilde = ( 2.0 * cell_dc * neig_dc ) /
                   ( neig_hxyz[xyz_idx] * cell_dc +
                     cell_hxyz[xyz_idx] * neig_dc );
          dhat = (J - dtilde * (phi_cell - phi_neig)) / (phi_neig + phi_cell);

          // Ensure positive definiteness
          if ((abs(dhat) > abs(dtilde)) && d_correct)
          {
            corrected = true;
            dtilde = 0.5 * J / phi_cell;
            dhat   = dtilde;
          }

          // Relax
          d_d_hat[g][cell][leak] = alpha * dhat
                                 + (1.0 - alpha) * d_d_hat[g][cell][leak];

          // Compute and set the off-diagonal matrix value.
          double val = (d_d_hat[g][cell][leak] - dtilde) / cell_hxyz[xyz_idx];

          int neig_row = neig_cell + g * d_group_size;
          flag = insert(row, neig_row, val, INSERT);
          Assert(flag);
        }

        // Compute leakage coefficient for this cell and surface.
        jo[leak] = dtilde + d_d_hat[g][cell][leak];

        //printf("Dhat(%2i) = %12.8f  Dtilde(%2i) = %12.8f \n", dhat, dtilde);

      } // leak loop

      // Net leakage coefficient.
      double jnet = (jo[1] + jo[0]) / d_mesh->dx(i) +
                    (jo[3] + jo[2]) / d_mesh->dy(j) +
                    (jo[5] + jo[4]) / d_mesh->dz(k);

      // Compute and set the diagonal matrix value.
      double val = jnet + cell_sr;

      flag = insert(row, row, val, INSERT);
      Assert(flag);

      // Add down/up scatter components
      {
        groups_t g_s =
          range<size_t>(d_material->lower(g, d_adjoint),
                        d_material->upper(g, d_adjoint),
                        true);
        groups_iter g_s_it = g_s.begin();
        for (; g_s_it != g_s.end(); ++g_s_it)
        {
          int gp = *g_s_it;

          // skip diagonal
          if (gp == g) continue;
          int col = cell + gp * d_group_size;

          double val = d_adjoint ? -d_material->sigma_s(m, gp, g)
                                 : -d_material->sigma_s(m, g, gp);
          flag = insert(row, col, val, INSERT);
          Assert(flag);
        }
      }
    } // row loop
  } // group loop
  if (corrected)
    printf("corrected!!!!!!!!!!!!!!!!\n");
  if (d_include_fission)
  {
    // Loop over all groups
    for (groups_iter g_it = d_groups.begin(); g_it != d_groups.end(); ++g_it)
    {
      int g  = *g_it;

      // Loop over all cells.
      for (int cell = 0; cell < d_group_size; cell++)
      {
        // Compute row index.
        int row = cell + g * d_group_size;

        // Define the data for this cell.
        int m = mat_map[cell];

        // Loop through source group.
        groups_iter gp_it = d_groups.begin();
        for (d_groups.begin(); gp_it != d_groups.end(); ++gp_it)
        {
          int gp  = *gp_it;

          // Compute column index.
          int col = cell + gp * d_group_size;

          // Fold the fission density with the spectrum.  Note that
          // we scale by keff and take the negative, since it's on the
          // left hand side.
          double v = 0.0;
          if (d_adjoint)
          {
            v = -d_material->nu_sigma_f(m, g) * d_material->chi(m, gp) / d_keff;
          }
          else
          {
            v = -d_material->nu_sigma_f(m, gp) * d_material->chi(m, g) / d_keff;
          }
          // Set the value. Note, we now have to add the value, since
          // in general fission contributes to nonzero cells.
          bool flag = insert(row, col, v, ADD);
          Assert(flag);
        }
      } // row loop
    } // group loop
  } // end fission block

  // Assemble.
  assemble();
}

//----------------------------------------------------------------------------//
template <class D>
double CMFDLossOperator<D>::albedo(const size_t side, const size_t g) const
{
  Require(side < 6);
  Require(g < d_number_groups);
  return d_albedo[side][g];
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class CMFDLossOperator<_1D>;
template class CMFDLossOperator<_2D>;
template class CMFDLossOperator<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file CMFDLossOperator.cc
//----------------------------------------------------------------------------//
