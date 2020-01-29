//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  WGDiffusionLossOperator.cc
 *  @brief WGDiffusionLossOperator member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "WGDiffusionLossOperator.hh"
#include <string>
#include <cmath>
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
WGDiffusionLossOperator::WGDiffusionLossOperator(SP_input    input,
                                                 SP_material material,
                                                 SP_mesh     mesh,
                                                 size_t      group)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_group(group)
  , d_albedo(6, 1.0)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_group < d_material->number_groups());

  // Set the dimension and group count
  d_dimension = d_mesh->dimension();

  // Set matrix dimensions
  Base::set_size(d_mesh->number_cells(), d_mesh->number_cells());

  // Nonzeros.  We have up to (diagonal + 2*dim neighbors) cells
  vec_int nnz(d_m, 1 + 2 * d_dimension);

  // Preallocate the matrix.
  preallocate(&nnz[0]);

  // Set the albedo based only on boundary condition.
  std::vector<std::string> boundary_name(6, "");
  boundary_name[Mesh::WEST]   = "bc_west";
  boundary_name[Mesh::EAST]   = "bc_east";
  boundary_name[Mesh::SOUTH]  = "bc_south";
  boundary_name[Mesh::NORTH]  = "bc_north";
  boundary_name[Mesh::BOTTOM] = "bc_bottom";
  boundary_name[Mesh::TOP]    = "bc_top";
  for (int b = 0; b < d_mesh->dimension() * 2; b++)
  {
    // Set the default to zero if it is an active boundary.  For 1/2-D problems,
    // we leave the "infinite" boundaries as reflective.
    d_albedo[b] = 0.0;
    if (d_input->check(boundary_name[b]))
      if (d_input->get<std::string>(boundary_name[b]) == "reflect")
        d_albedo[b] = 1.0;
  }

  // Build the matrix.
  build();

  //print_matlab("WGDSA.out");

}

//---------------------------------------------------------------------------//
void WGDiffusionLossOperator::construct()
{
  build();
}

//---------------------------------------------------------------------------//
void WGDiffusionLossOperator::build()
{
  using std::cout;
  using std::endl;

  // The matrix dimension is the total number of cells.
  size_t size = d_mesh->number_cells();

  // Get the material map.
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  // Error flag
  bool flag;

  double tmp = 1.0;

  // Loop over all matrix rows, which, because of the ordering,
  // is the same as the cell index.
  for (int row = 0; row < size; row++)
  {

    // Define the data for this cell.
    int m = mat_map[row];

    double cell_dc = tmp * d_material->diff_coef(m, d_group);
    Assert(cell_dc > 0.0);


    double cell_sr = d_material->sigma_t(m, d_group) -
                     d_material->sigma_s(m, d_group, d_group);
    //cell_sr *= 2.0;
//    double cell_sr = d_material->sigma_a(m, d_group);
//    for (size_t g = 0; g < d_material->number_groups(); ++g)
//    {
//      if (g == d_group) continue;
//      cell_sr += d_material->sigma_s(m, g, d_group);
//    }

    // Get the directional indices.
    int i = d_mesh->cell_to_i(row);
    int j = d_mesh->cell_to_j(row);
    int k = d_mesh->cell_to_k(row);

    // Cell width vector.
    double cell_hxyz[3] = {d_mesh->dx(i), d_mesh->dy(j), d_mesh->dz(k)};

    // Direction-specific leakage coefficients.
    double jo[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Index arrays to help determine if a cell surface is on the boundary.
    int bound[6] = {i, i, j, j, k, k};
    int nxyz[3][2] = {{0, d_mesh->number_cells_x()-1},
                      {0, d_mesh->number_cells_y()-1},
                      {0, d_mesh->number_cells_z()-1}};

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

        dtilde = ( 2.0 * cell_dc * (1.0 - d_albedo[leak]) ) /
                 ( 4.0 * cell_dc * (1.0 + d_albedo[leak]) +
                   (1.0 - d_albedo[leak]) * cell_hxyz[xyz_idx] );

      }
      else // not a boundary
      {

        // Get the neighbor data.
        int neig_row = d_mesh->index(neig_idx[0], neig_idx[1], neig_idx[2]);
        Assert(neig_row >= 0);
        int ii = d_mesh->cell_to_i(neig_row);
        int jj = d_mesh->cell_to_j(neig_row);
        int kk = d_mesh->cell_to_k(neig_row);

        // Neighbor volume and diffusion coefficient.
        double neig_hxyz[3] = {d_mesh->dx(ii), d_mesh->dy(jj), d_mesh->dz(kk)};
        double neig_dc = tmp * d_material->diff_coef(mat_map[neig_row], d_group);

        // Compute dtilde.
        dtilde = ( 2.0 * cell_dc * neig_dc ) /
                 ( neig_hxyz[xyz_idx] * cell_dc +
                   cell_hxyz[xyz_idx] * neig_dc );

        // Compute and set the off-diagonal matrix value.
        double val = - dtilde / cell_hxyz[xyz_idx];
        flag = insert(row, neig_row, val, INSERT);

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
   flag = insert(row, row, val, INSERT);

  } // row loop

  // Assemble.
  assemble();

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file WGDiffusionLossOperator.cc
//---------------------------------------------------------------------------//
