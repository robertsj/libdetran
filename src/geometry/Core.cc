/*
 * Core.cc
 *
 *  Created on: Apr 16, 2012
 *      Author: robertsj
 */

// Detran
#include "Core.hh"

// System
#include <cmath>
#include <iostream>

namespace detran
{

Core::Core(int dimension, vec_assembly assemblies, vec_int assembly_map)
  : d_dimension(dimension)
  , d_assemblies(assemblies)
{
  Require(d_dimension > 0);
  Require(d_assemblies.size() > 0);
  finalize(assembly_map);
}

Core::Core(int dimension)
 : d_dimension(dimension)
{
  Require(d_dimension > 0);
}

void Core::add_assembly(SP_assembly assembly)
{
  Require(assembly);
  d_assemblies.push_back(assembly);
}

void Core::finalize(vec_int assembly_map)
{
  using std::cout;
  using std::endl;

  Require(assembly_map.size() == d_dimension*d_dimension);
  d_assembly_map = assembly_map;

  // Set number of cells.  This *assumes* all pins have the same meshing.
  int number_row = std::sqrt(assembly_map.size());
  Assert(assembly_map.size() == number_row*number_row);
  int number_cells_x = (d_assemblies[0]->mesh())->number_cells_x() * number_row;
  int number_cells_y = (d_assemblies[0]->mesh())->number_cells_x() * number_row;
  int number_cells = number_cells_x * number_cells_y;

  // Compute the widths.
  double width = d_assemblies[0]->mesh()->dx(0);

  // Fine mesh edges
  vec_dbl edges(number_cells_x + 1, 0.0);

  // Temporary maps.
  vec_int tmp_mat_map(number_cells, 0);
  vec_int tmp_reg_map(number_cells, 0);
  vec_int tmp_pin_map(number_cells, 0);

  int assembly_count = 1;
  int j_save = 0;
  int i_save = 0;

  for (int j = 0; j < number_row; j++)
  {
    // Fine mesh y range for this jth coarse mesh.
    int j1 = j_save;
    int j2 = j_save + d_assemblies[0]->mesh()->number_cells_y();

    // Do edges once, since they are shared.
    for (int jj = j1; jj < j2; jj++)
    {
      // All pins have same mesh in a direction.
      edges[jj + 1] = edges[jj] + d_assemblies[assembly_map[j*number_row]]->mesh()->dy(j1 - j_save);
    }

    for (int i = 0; i < number_row; i++)
    {

      // Pin index
      int pin = i + j*number_row;

      // Fine mesh x range.
      int i1 = i_save;
      int i2 = i_save + d_assemblies[0]->mesh()->number_cells_x();

      // Get the material and region maps for this pin.
      vec_int m = d_assemblies[assembly_map[pin]]->mesh()->mesh_map("MATERIAL");
      vec_int r = d_assemblies[assembly_map[pin]]->mesh()->mesh_map("REGION");

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
          tmp_pin_map[cell] = assembly_count;
          count++;
        }
      }
      i_save += d_assemblies[0]->mesh()->number_cells_x();
      assembly_count++;
    }
    i_save = 0;
    j_save += d_assemblies[0]->mesh()->number_cells_y();

  }

  // Create my mesh.
  d_mesh = new Mesh2D(edges, edges, tmp_mat_map);
  // Add maps.
  d_mesh->add_mesh_map("REGION", tmp_reg_map);
  // Assigns unique edit region for each pin in an Core.
  d_mesh->add_mesh_map("ASSEMBLIES", tmp_pin_map);

}

} // end namespace detran







