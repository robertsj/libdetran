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
#include <iostream>

namespace detran
{

Assembly::Assembly(int dimension, vec_pincell pincells, vec_int pincell_map)
  : d_dimension(dimension)
  , d_pincells(pincells)
  , d_pincell_map(pincell_map)
{
  Require(d_dimension > 0);
  Require(d_pincells.size() > 0);
  finalize(d_pincell_map);
}

Assembly::Assembly(int dimension)
 : d_dimension(dimension)
{
  Require(d_dimension > 0);
}

void Assembly::add_pincell(SP_pincell pin)
{
  Require(pin);
  d_pincells.push_back(pin);
}

void Assembly::finalize(vec_int pincell_map)
{
  using std::cout;
  using std::endl;

  // Verify that pin cells are consistent, using first
  // pin as the reference.
  for (int i = 1; i < d_pincells.size(); i++)
  {
    Insist(d_pincells[0]->mesh()->number_cells_x() ==
           d_pincells[i]->mesh()->number_cells_x(),
           "X fine mesh inconsistent.");
    Insist(d_pincells[0]->mesh()->number_cells_y() ==
           d_pincells[i]->mesh()->number_cells_y(),
           "Y fine mesh inconsistent.");
  }

  Require(pincell_map.size() == d_dimension*d_dimension);
  d_pincell_map = pincell_map;

  // Set number of cells.  This *assumes* all pins have the same meshing.
  int number_pins_row = std::sqrt(pincell_map.size());
  Assert(pincell_map.size() == number_pins_row*number_pins_row);
  int number_cells_x =
    (d_pincells[0]->mesh())->number_cells_x() * number_pins_row;
  int number_cells_y =
    (d_pincells[0]->mesh())->number_cells_x() * number_pins_row;
  int number_cells = number_cells_x * number_cells_y;

  // Compute the widths.
  double width = d_pincells[0]->mesh()->dx(0);

  // Fine mesh edges
  vec_dbl edges(number_cells_x + 1, 0.0);

  // Temporary maps.
  vec_int tmp_mat_map(number_cells, 0);
  vec_int tmp_reg_map(number_cells, 0);
  vec_int tmp_pin_map(number_cells, 0);

  int pin_count = 1;
  int j_save = 0;
  int i_save = 0;

  for (int j = 0; j < number_pins_row; j++)
  {
    // Fine mesh y range for this jth coarse mesh.
    int j1 = j_save;
    int j2 = j_save + d_pincells[0]->mesh()->number_cells_y();

    // Do edges once, since they are shared.
    for (int jj = j1; jj < j2; jj++)
    {
      // All pins have same mesh in a direction.
      edges[jj + 1] = edges[jj] +
        d_pincells[pincell_map[j*number_pins_row]]->mesh()->dy(jj - j_save);
    }

    for (int i = 0; i < number_pins_row; i++)
    {

      // Pin index
      int pin = i + j*number_pins_row;

      // Fine mesh x range.
      int i1 = i_save;
      int i2 = i_save + d_pincells[0]->mesh()->number_cells_x();

      // Get the material and region maps for this pin.
      vec_int m = d_pincells[pincell_map[pin]]->mesh()->mesh_map("MATERIAL");
      vec_int r = d_pincells[pincell_map[pin]]->mesh()->mesh_map("REGION");

      // Assign the values.
      int count = 0;
      for (int jj = j1; jj < j2; jj++)
      {
        //cout << " jj = " << jj << endl;
        for (int ii = i1; ii < i2; ii++)
        {
          int cell = ii + jj*number_cells_x;
          tmp_mat_map[cell] = m[count];
          tmp_reg_map[cell] = r[count];
          tmp_pin_map[cell] = pin_count - 1;
          count++;
        }
      }
      i_save += d_pincells[0]->mesh()->number_cells_x();
      pin_count++;
    }
    i_save = 0;
    j_save += d_pincells[0]->mesh()->number_cells_y();

  }

  // Create my mesh.
  d_mesh = new Mesh2D(edges, edges, tmp_mat_map);
  // Add maps.
  d_mesh->add_mesh_map("REGION", tmp_reg_map);
  // Assigns unique edit region for each pin in an assembly.
  d_mesh->add_mesh_map("PINS", tmp_pin_map);

}

} // end namespace detran



