/*
 * Assembly.cc
 *
 *  Created on: Apr 14, 2012
 *      Author: robertsj
 */

// Detran
#include "Assembly.hh"

// System
#include <cmath>

namespace detran
{

Assembly::Assembly(vec_pincell pincells, vec_int pincell_map)
  : d_pincells(pincells)
  , d_pincell_map(pincell_map)
{

  // Set number of cells.  This *assumes* all pins have the same
  // meshing, though not necessarily uniform.
  int number_pins_row = std::sqrt(pincell_map.size());
  d_number_cells_x = d_pincells[0]->number_cells_x() * number_pins_row;
  d_number_cells_y = d_pincells[0]->number_cells_x() * number_pins_row;
  d_number_cells_z = 1;
  d_number_cells = d_number_cells_x * d_number_cells_y;

  // Compute the widths.
  double width = d_pincells[0]->dx(0);
  d_dx = vec_dbl(d_number_cells_x, width);
  d_dy = vec_dbl(d_number_cells_y, width);

  vec_int tmp_mat_map(d_number_cells, 0);
  vec_int tmp_reg_map(d_number_cells, 0);
  vec_int tmp_pin_map(d_number_cells, 0);

  int pin_count = 0;
  int j_save = 0;
  int i_save = 0;
  for (int j = 0; j < number_pins_row; j++)
  {
    // Fine mesh y range for this jth coarse mesh.
    int j1 = j_save;
    int j2 = j_save + d_pincells[0]->number_cells_y();

    for (int jj = j1; jj < j2; jj++)
    {
      // All pins have same mesh in a direction.
      d_dy[jj] = d_pincells[j*number_pins_row]->dy(j1 - j_save);
    }

    for (int i = 0; i < number_pins_row; i++)
    {

      // Pin index
      int pin = i + j*number_pins_row;

      // Fine mesh x range.
      int i1 = i_save;
      int i2 = i_save + d_pincells[0]->number_cells_x();

      // Set fine mesh width in x direction.
      if (j == 0)
      {
        for (int ii = i1; ii < i2; ii++)
        {
          // All pins have same mesh in a row.
          d_dx[ii] = d_pincells[i]->dx(i1 - i_save);
        }
      }

      // Get the material and region maps for this pin.
      vec_int m = d_pincells[pin]->mesh_map("MATERIAL");
      vec_int r = d_pincells[pin]->mesh_map("REGION");

      // Assign the values.
      int count = 0;
      for (int jj = j1; jj < j2; jj++)
      {
        for (int ii = i1; ii < i2; ii++)
        {
          tmp_mat_map[index(ii, jj, 0)] = m[count];
          tmp_reg_map[index(ii, jj, 0)] = r[count];
          tmp_pin_map[index(ii, jj, 0)] = pin_count;
          count++;
        }
      }
      pin_count++;
    }

  }

  // Add maps.
  add_mesh_map("MATERIAL", tmp_mat_map);
  add_mesh_map("REGION",   tmp_reg_map);
  // Assigns unique edit region for each pin in an assembly.
  add_mesh_map("PINS",     tmp_pin_map);

}


} // end namespace detran



