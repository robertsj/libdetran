//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DiffusionLossOperator.cc
 *  @brief  DiffusionLossOperator
 *  @author Jeremy Roberts
 *  @date   Sep 10, 2012
 */
//---------------------------------------------------------------------------//

#include "DiffusionLossOperator.hh"
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
DiffusionLossOperator::DiffusionLossOperator(SP_input       input,
                                             SP_material    material,
                                             SP_mesh        mesh,
                                             const bool     include_fission,
                                             const size_t   cutoff,
                                             const bool     adjoint,
                                             const double   keff)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_include_fission(include_fission)
  , d_group_cutoff(cutoff)
  , d_adjoint(adjoint)
  , d_keff(keff)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Set the dimension and group count
  d_dimension = d_mesh->dimension();
  d_number_groups = d_material->number_groups();
  d_number_active_groups = d_number_groups - d_group_cutoff;
  d_group_size    = d_mesh->number_cells();

  // Set matrix dimensions
  Base::set_size(d_number_active_groups*d_group_size,
                 d_number_active_groups*d_group_size);

  // Default albedos to 1.0.  For dimensions in play, this will be
  // overwritten by the default boundary.
  d_albedo.resize(6,  vec_dbl(d_number_active_groups, 1.0));

  // Nonzeros.  We have up to
  //   diagonal + 2*dim neighbors + num_groups coupling from scatter/fission
  vec_int nnz(d_m, 1 + 2 * d_dimension + d_number_active_groups);

  // Preallocate the matrix.  Note, PETSc documentation suggests getting
  // this right is extremely important.
  preallocate(&nnz[0]);

  // Set the albedo.  First, check if the input has an albedo
  // entry.  If it does, this is the default way to set the
  // condition.  Otherwise, check the boundary conditions.
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
    for (int g = d_group_cutoff; g < d_number_groups; g++)
    {
      for (int b = 0; b < d_mesh->dimension() * 2; b++)
      {
        d_albedo[b][g - d_group_cutoff] = 0.0;
        if (d_input->check(boundary_name[b]))
          if (d_input->get<std::string>(boundary_name[b]) == "reflect")
            d_albedo[b][g - d_group_cutoff] = 1.0;
      }
    }
  }

  // Build the matrix with the initial keff guess.
  build();

}

//---------------------------------------------------------------------------//
void DiffusionLossOperator::construct(const double keff)
{
  d_keff = keff;
  build();
}



//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void DiffusionLossOperator::build()
{

  using std::cout;
  using std::endl;

  // Get the material map.
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  d_material->display();

  bool db = false;

  for (int gg = 0; gg < d_number_active_groups; gg++)
  {
    int g = gg + d_group_cutoff;

    // Loop over all cells.
    for (int cell = 0; cell < d_group_size; cell++)
    {
      if (db) cout << "  cell = " << cell << endl;

      // Compute row index.
      int row = cell + gg * d_group_size;

      if (db) cout << "    row = " << row << endl;

      // Define the data for this cell.
      size_t m = mat_map[cell];

      double cell_dc = d_material->diff_coef(m, g);
      std::cout << " m=" << m << " g=" << g << " dc=" << cell_dc << std::endl;
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


      // Flag for insertion status.  If false, it means we didn't add
      // the value and should counter accordingly (here, we assert)
      bool flag;

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

        // \todo It might make sense to template the coupling coefficient

        // Compute coupling coefficient
        double dtilde = 0.0;
        if (bound[leak] == nxyz[xyz_idx][dir_idx])
        {

          dtilde = ( 2.0 * cell_dc * (1.0 - d_albedo[leak][gg]) ) /
                   ( 4.0 * cell_dc * (1.0 + d_albedo[leak][gg]) +
                    (1.0 - d_albedo[leak][gg]) * cell_hxyz[xyz_idx] );

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
          int neig_row = neig_cell + gg * d_group_size;

          if (db) cout << "      col = " << neig_row << " " << neig_cell << endl;

          flag = insert(row, neig_row, val, INSERT);
          Assert(flag);
        }

        // Compute leakage coefficient for this cell and surface.
        jo[leak] = dtilde;

      } // leak loop

      // Net leakage coefficient.
      double jnet = (jo[1] + jo[0]) / d_mesh->dx(i) +
                    (jo[3] + jo[2]) / d_mesh->dy(j) +
                    (jo[5] + jo[4]) / d_mesh->dz(k);

     // Compute and set the diagonal matrix value.
     double val = jnet + cell_sr;
     //val = 1.0;
     if (db) cout << "      col = " << row << " " << gg <<  endl;
    // cout << " jnet=" << jnet << " cellsr=" << cell_sr << endl;
     flag = insert(row, row, val, INSERT);
     Assert(flag);

     // Add downscatter component.
     int lower = d_material->lower(g);
     if (d_group_cutoff > lower) lower = d_group_cutoff;
     for (int gp = lower; gp < g; gp++)
     {
       int col = cell + (gp - d_group_cutoff) * d_group_size;
       if (db) cout << "  ds  col = " << col << endl;
       double val = -d_material->sigma_s(m, g, gp);
       if (!d_adjoint)
         flag = insert(row, col, val, INSERT);
       else
         flag = insert(col, row, val, INSERT);
       Assert(flag);
     }

     // Add upscatter component.
     for (int gp = g + 1; gp <= d_material->upper(g); gp++)
     {
       int col = cell + (gp - d_group_cutoff) * d_group_size;
       double val = -d_material->sigma_s(m, g, gp);
       if (db) cout << "  us  col = " << col << endl;
       if (!d_adjoint)
         flag = insert(row, col, val, INSERT);
       else
         flag = insert(col, row, val, INSERT);
       Assert(flag);
     }

    } // row loop

  } // group loop

  if (d_include_fission)
  {
    // Loop over all groups
    for (int g = d_group_cutoff; g < d_number_groups; g++)
    {
      // Loop over all cells.
      for (int cell = 0; cell < d_group_size; cell++)
      {
        // Compute row index.
        int row = cell + (g - d_group_cutoff) * d_group_size;

        // Define the data for this cell.
        int m = mat_map[cell];

        // Get the directional indices.
        int i = d_mesh->cell_to_i(cell);
        int j = d_mesh->cell_to_j(cell);
        int k = d_mesh->cell_to_k(cell);

        // Loop through source group.
        for (int gp = d_group_cutoff; gp < d_number_groups; gp++)
        {
          // Compute column index.
          int col = cell + (gp - d_group_cutoff) * d_group_size;
          if (db) cout << "      col = " << col << endl;
          // Fold the fission density with the spectrum.  Note that
          // we scale by keff and take the negative, since it's on the
          // left hand side.
          double val = -d_material->nu_sigma_f(m, gp) *
                        d_material->chi(m, g) / d_keff;

          // Set the value. Note, we now have to add the value, since
          // in general fission contributes to nonzero cells.
          bool flag;
          if (!d_adjoint)
            flag = insert(row, col, val, ADD);
          else
            flag = insert(col, row, val, ADD);
          Assert(flag);
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
