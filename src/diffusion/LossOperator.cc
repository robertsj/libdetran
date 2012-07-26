//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LossOperator.cc
 * \brief  LossOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC

// Detran
#include "LossOperator.hh"

// System
#include <string>
#include <cmath>
#include <iostream>

namespace detran_diffusion
{

LossOperator::LossOperator(SP_input    input,
                           SP_material material,
                           SP_mesh     mesh)
  : BaseOperator(input, material, mesh)
  , d_albedo(6, vec_dbl(d_material->number_groups(), 1.0))
  , d_number_groups(d_material->number_groups())
  , d_group_size(d_mesh->number_cells())
{
  // Preconditions
  using std::string;

  // The matrix dimension
  d_size = d_group_size * d_number_groups;


  // Set the albedo.  First, check if the input has an albedo
  // entry.  If it does, this is the default way to set the
  // condition.  Otherwise, check the boundary conditions.
  if (d_input->check("albedo"))
  {
//    vec_dbl albedo = d_input->get<vec_dbl>("albedo");
//    Insist(albedo.size() == d_albedo.size(),
//      "The user-defined albedo vector is the wrong size.");
//    d_albedo = albedo;
  }
  else
  {
    for (int g = 0; g < d_material->number_groups(); g++)
    {
      // Left and right
      d_albedo[detran::Mesh::LEFT][g] = 0.0;
      if (d_input->check("bc_left"))
        if (d_input->get<string>("bc_left") == "reflect")
          d_albedo[detran::Mesh::LEFT][g] = 1.0;
      d_albedo[detran::Mesh::RIGHT][g] = 0.0;
      if (d_input->check("bc_right"))
        if (d_input->get<string>("bc_right") == "reflect")
          d_albedo[detran::Mesh::RIGHT][g] = 1.0;

      if (d_dimension > 1)
      {
        // Bottom and top
        d_albedo[detran::Mesh::BOTTOM][g] = 0.0;
        if (d_input->check("bc_bottom"))
          if (d_input->get<string>("bc_bottom") == "reflect")
            d_albedo[detran::Mesh::BOTTOM][g] = 1.0;
        d_albedo[detran::Mesh::TOP][g] = 0.0;
        if (d_input->check("bc_top"))
          if (d_input->get<string>("bc_top") == "reflect")
            d_albedo[detran::Mesh::TOP][g] = 1.0;
      }

      if (d_dimension == 3)
      {
        // South and north
        d_albedo[detran::Mesh::SOUTH][g] = 0.0;
        if (d_input->check("bc_south"))
          if (d_input->get<string>("bc_south") == "reflect")
            d_albedo[detran::Mesh::SOUTH][g] = 1.0;
        d_albedo[detran::Mesh::NORTH][g] = 0.0;
        if (d_input->check("bc_north"))
          if (d_input->get<string>("bc_north") == "reflect")
            d_albedo[detran::Mesh::NORTH][g] = 1.0;
      }
    }
  }

  // Construct the matrix.
  construct();

}


//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

void LossOperator::construct()
{

  using std::cout;
  using std::endl;

  // Create the matrix.
  MatCreate(PETSC_COMM_SELF, &d_operator);
  MatSetType(d_operator, MATAIJ);
  MatSetSizes(d_operator, PETSC_DECIDE, PETSC_DECIDE, d_size, d_size);

  // Get the material map.
  detran::vec_int mat_map = d_mesh->mesh_map("MATERIAL");

  for (int g = 0; g < d_number_groups; g++)
  {

    //cout << " g = " << g << endl;

    // Loop over all cells.
    for (int cell = 0; cell < d_group_size; cell++)
    {
      //cout << " cell = " << cell << endl;

      // Compute row index.
      int row = cell + g * d_group_size;

      //cout << " row = " << row << endl;

      // Define the data for this cell.
      int m = mat_map[cell];

      //cout << " m = " << m << endl;

      double cell_dc = d_material->diff_coef(m, g);
      Assert(cell_dc > 0.0);
      double cell_sr = d_material->sigma_t(m, g) -
                       d_material->sigma_s(m, g, g);

      // Get the directional indices.
      int i = d_mesh->cell_to_i(cell);
      int j = d_mesh->cell_to_j(cell);
      int k = d_mesh->cell_to_k(cell);

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
          int neig_cell = d_mesh->index(neig_idx[0], neig_idx[1], neig_idx[2]);
          Assert(neig_cell >= 0);
          int ii = d_mesh->cell_to_i(neig_cell);
          int jj = d_mesh->cell_to_j(neig_cell);
          int kk = d_mesh->cell_to_k(neig_cell);

          //cout << " neig_row " << neig_row << endl;

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

     // Add downscatter component.
     for (int gp = d_material->lower(g); gp < g; gp++)
     {
       int col = cell + gp * d_group_size;
       double val = -d_material->sigma_s(m, g, gp);
       MatSetValue(d_operator, row, col, val, INSERT_VALUES);
     }

     // Add upscatter component.
     for (int gp = g + 1; gp <= d_material->upper(g); gp++)
     {
       int col = cell + gp * d_group_size;
       double val = -d_material->sigma_s(m, g, gp);
       MatSetValue(d_operator, row, col, val, INSERT_VALUES);
     }

     // Add the fission component, if this is a fixed source
     // multiplying problem.
     if (1==0)
     {
       // finish me
     }

    } // row loop

  } // group loop

  // Assemble.
  MatAssemblyBegin(d_operator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(d_operator, MAT_FINAL_ASSEMBLY);

}

} // end namespace detran_diffusion

#endif

//---------------------------------------------------------------------------//
//              end of file LossOperator.cc
//---------------------------------------------------------------------------//
