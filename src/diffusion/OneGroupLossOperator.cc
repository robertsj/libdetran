//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OneGroupLossOperator.cc
 * \brief  OneGroupLossOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.hh"

#ifdef DETRAN_ENABLE_PETSC

// Detran
#include "OneGroupLossOperator.hh"

// System
#include <string>
#include <cmath>
#include <iostream>

namespace detran_diffusion
{

OneGroupLossOperator::OneGroupLossOperator(SP_input    input,
                                           SP_material material,
                                           SP_mesh     mesh,
                                           int         group)
  : BaseOperator(input, material, mesh)
  , d_albedo(6, 1.0)
  , d_group(group)
{
  // Preconditions
  Require(d_group >= 0);
  Require(d_group < d_material->number_groups());

  using std::string;

  // Set the albedo.  First, check if the input has an albedo
  // entry.  If it does, this is the default way to set the
  // condition.  Otherwise, check the boundary conditions.
  if (d_input->check("albedo"))
  {
    vec_dbl albedo = d_input->get<vec_dbl>("albedo");
    Insist(albedo.size() == d_albedo.size(),
      "The user-defined albedo vector is the wrong size.");
    d_albedo = albedo;
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
    for (int b = 0; b < d_mesh->dimension() * 2; b++)
    {
      // Set the default to zero if it is an active boundary.  For 1/2-D problems,
      // we leave the "infinite" boundaries as reflective.
      d_albedo[b] = 0.0;
      if (d_input->check(boundary_name[b]))
        if (d_input->get<string>(boundary_name[b]) == "reflect")
          d_albedo[b] = 1.0;
    }
  }

  // Construct the matrix.
  construct();

}


// IMPLEMENTATION

//
void OneGroupLossOperator::construct()
{
  using std::cout;
  using std::endl;

  // The matrix dimension is the total number of cells.
  d_size = d_mesh->number_cells();

  // Create the matrix.
  MatCreate(PETSC_COMM_SELF, &d_operator);
  MatSetType(d_operator, MATAIJ);
  MatSetSizes(d_operator, PETSC_DECIDE, PETSC_DECIDE, d_size, d_size);

  // Get the material map.
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  // Loop over all matrix rows, which, because of the ordering,
  // is the same as the cell index.
  for (int row = 0; row < d_size; row++)
  {
    //cout << " row = " << row << endl;

    // Define the data for this cell.
    int m = mat_map[row];
    //cout << " m = " << m << " g = " << d_group << endl;

    double cell_dc = d_material->diff_coef(m, d_group);
    Assert(cell_dc > 0.0);
    double cell_sr = d_material->sigma_t(m, d_group) -
                     d_material->sigma_s(m, d_group, d_group);

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
        double neig_dc = d_material->diff_coef(mat_map[neig_row], d_group);

        // Compute dtilde.
        dtilde = ( 2.0 * cell_dc * neig_dc ) /
                 ( neig_hxyz[xyz_idx] * cell_dc +
                   cell_hxyz[xyz_idx] * neig_dc );

        // Compute and set the off-diagonal matrix value.
        double val = - dtilde / cell_hxyz[xyz_idx];
        MatSetValue(d_operator, row, neig_row, val, INSERT_VALUES);

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
   MatSetValue(d_operator, row, row, val, INSERT_VALUES);

  } // row loop

  // Assemble.
  MatAssemblyBegin(d_operator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(d_operator, MAT_FINAL_ASSEMBLY);

}

} // end namespace detran

#endif

//---------------------------------------------------------------------------//
//              end of file OneGroupLossOperator.cc
//---------------------------------------------------------------------------//
