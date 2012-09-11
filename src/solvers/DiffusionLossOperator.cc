//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionLossOperator.cc
 * \brief  DiffusionLossOperator 
 * \author Jeremy Roberts
 * \date   Sep 10, 2012
 */
//---------------------------------------------------------------------------//

#include "DiffusionLossOperator.hh"
#include <iostream>

namespace detran
{

DiffusionLossOperator::DiffusionLossOperator(SP_input       input,
                                             SP_material    material,
                                             SP_mesh        mesh,
                                             const bool     include_fission,
                                             const double   keff)
  : OperatorMatrix(material->number_groups()*mesh->number_cells(),
                   material->number_groups()*mesh->number_cells())
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_dimension(mesh->dimension())
  , d_include_fission(false)
  , d_keff(keff)
{

  // Nonzeros.  We have up to the number of groups in a block row
  // due to scattering and potentially fission.  Additionally,
  // we have 2*dimension entries off the energy block due to neighbors.
  vec_int nnz(d_number_rows, 2 * d_dimension + d_number_groups);

  // Preallocate the matrix.  Note, PETSc documentation suggests getting
  // this right is extremely important.
  preallocate(nnz);

  // Build the matrix with the initial keff guess.
  build();

}

void DiffusionLossOperator::construct(const double keff)
{
  d_keff = keff;
  build();
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

/*!
 *  \page mesh-centered-diffusion Mesh-Centered Finite-Difference Diffusion
 *
 *  In this note, we describe our implementation of finite difference
 *  diffusion.  We assume a Cartesian grid with cells in which the materials,
 *  flux, and source are taken to be constant.
 *
 *
 *
 */
void DiffusionLossOperator::build()
{

  using std::cout;
  using std::endl;

  // Get the material map.
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  for (int g = 0; g < d_number_groups; g++)
  {
    // Loop over all cells.
    for (int cell = 0; cell < d_group_size; cell++)
    {
      // Compute row index.
      int row = cell + g * d_group_size;

      // Define the data for this cell.
      size_t m = mat_map[cell];

      double cell_dc = d_material->diff_coef(m, g);
      Assert(cell_dc > 0.0);
      double cell_sr = d_material->sigma_t(m, g) -
                       d_material->sigma_s(m, g, g);

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
      int nxyz[3][2] = {0, d_mesh->number_cells_x()-1,
                        0, d_mesh->number_cells_y()-1,
                        0, d_mesh->number_cells_z()-1};


      // For each spatial cell, there are 6 faces that connect the
      // cell to a neighbor or the global boundary.  Looping through
      // each of the faces, we can use the indices of the cell to
      // determine whether the face is between cells or is a boundary
      // surface.

      // leak --> 0=-x, 1=+x, 2=-y, 3=+y, 4=-z, 5=+z
      for (int leak = 0; leak < 6; leak++)
      {

        // Determine whether this is a left/right, bottom/top,
        // or south/north boundary.
        int xyz_idx = std::floor(leak / 2);

        // Determine the direction, e.g. left (-) vs right (+).
        int dir_idx = 1 - ((leak + 1) % 2);

        // Put i, j, k into an array.  The neighbor indices are
        // a perturbation of this.
        int neig_idx[3] = {i, j, k};

        // Determine whether the neighbor is positive (+1) or negative (-1)
        // relative to the surface under consideration, and then decrement
        // the appropriate x, y, or z index.
        int shift_idx   = -2 * ((leak + 1) % 2) + 1;
        neig_idx[xyz_idx] += shift_idx;

        // Compute coupling coefficient
        double dtilde = 0.0;
        if (bound[leak] == nxyz[xyz_idx][dir_idx])
        {

          dtilde = ( 2.0 * cell_dc * (1.0 - d_albedo[leak][g]) ) /
                   ( 4.0 * cell_dc * (1.0 + d_albedo[leak][g]) +
                     (1.0 - d_albedo[leak][g]) * cell_hxyz[xyz_idx] );

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

          // Compute and set the off-diagonal matrix value.
          double val = - dtilde / cell_hxyz[xyz_idx];
          int neig_row = neig_cell + g * d_group_size;
          insert_values(1, &row, 1, &neig_row, &val, INSERT_VALUES);
        }

        // Compute leakage coefficient for this cell and surface.
        jo[leak] = (double)shift_idx * dtilde;

      } // leak loop

      // Net leakage coefficient.
      double jnet = (jo[1] - jo[0]) / d_mesh->dx(i) +
                    (jo[3] - jo[2]) / d_mesh->dy(j) +
                    (jo[5] - jo[4]) / d_mesh->dz(k);

     // Compute and set the diagonal matrix value.
     double val = jnet + cell_sr;
     insert_values(1, &row, 1, &row, &val, INSERT_VALUES);

     // Add downscatter component.
     for (int gp = d_material->lower(g); gp < g; gp++)
     {
       int col = cell + gp * d_group_size;
       double val = -d_material->sigma_s(m, g, gp);
       insert_values(1, &row, 1, &col, &val, INSERT_VALUES);
     }

     // Add upscatter component.
     for (int gp = g + 1; gp <= d_material->upper(g); gp++)
     {
       int col = cell + gp * d_group_size;
       double val = -d_material->sigma_s(m, g, gp);
       insert_values(1, &row, 1, &col, &val, INSERT_VALUES);
     }

    } // row loop

  } // group loop

  if (d_include_fission)
  {
    // Loop over all groups
    for (int g = 0; g < d_number_groups; g++)
    {
      // Loop over all cells.
      for (int cell = 0; cell < d_group_size; cell++)
      {
        // Compute row index.
        int row = cell + g * d_group_size;

        // Define the data for this cell.
        int m = mat_map[cell];

        // Get the directional indices.
        int i = d_mesh->cell_to_i(cell);
        int j = d_mesh->cell_to_j(cell);
        int k = d_mesh->cell_to_k(cell);

        // Loop through source group.
        for (int gp = 0; gp < d_number_groups; gp++)
        {
          // Compute column index.
          int col = cell + gp * d_group_size;

          // Fold the fission density with the spectrum.  Note that
          // we scale by keff and take the negative, since it's on the
          // left hand side.
          double val = -d_material->nu_sigma_f(m, gp) *
                       d_material->chi(m, g) / d_keff;

          // Set the value.
          insert_values(1, &row, 1, &col, &val, ADD_VALUES);
        }

      } // row loop

    } // group loop

  } // end fission block

  // Assemble.
  assemble();

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file DiffusionLossOperator.cc
//---------------------------------------------------------------------------//
