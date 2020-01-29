//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Core.cc
 *  @brief Core class member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "Core.hh"
#include <cmath>
#include <iostream>

namespace detran_geometry
{

//---------------------------------------------------------------------------//
Core::Core(int dimension, vec_assembly assemblies, vec_int assembly_map)
  : d_number_x(dimension)
  , d_number_y(dimension)
  , d_assemblies(assemblies)
{
  Require(d_number_x > 0);
  Require(d_assemblies.size() > 0);
  finalize(assembly_map);
}

//---------------------------------------------------------------------------//
Core::Core(int dimension)
 : d_number_x(dimension)
 , d_number_y(dimension)

{
  Require(d_number_x > 0);
}

//---------------------------------------------------------------------------//
void Core::add_assembly(SP_assembly assembly)
{
  Require(assembly);
  d_assemblies.push_back(assembly);
}

//---------------------------------------------------------------------------//
void Core::finalize(vec_int assembly_map)
{
  using std::cout;
  using std::endl;

  Insist(assembly_map.size() == d_number_x*d_number_y,
    "Assembly map is wrong size.");
  d_assembly_map = assembly_map;

  // Set number of cells.  This *assumes* all pins have the same meshing.
  int number_row = std::sqrt(assembly_map.size());
  Assert(assembly_map.size() == number_row*number_row);
  int number_cells_x =
    (d_assemblies[0]->mesh())->number_cells_x() * number_row;
  int number_cells_y =
    (d_assemblies[0]->mesh())->number_cells_x() * number_row;
  int number_cells = number_cells_x * number_cells_y;

  // Number of pins per dimension
  d_number_assemblies   = d_number_x*d_number_y;
  d_number_pincells     =
    d_number_assemblies * d_assemblies[0]->number_pincells();

  int ass_dim = d_assemblies[0]->dimension();


  // Create map from the assembly pin indices to the core level.
  vec_int core_pins(d_number_pincells, -1);
  // Loop through all pins.
  for (int j = 0; j < d_number_y; j++)
  {
    for (int i = 0; i < d_number_x; i++)
    {
      // What pin?
      int pin = i + j * d_number_x;

      // What assembly?
      int a_i = std::floor(i / ass_dim);
      int a_j = std::floor(j / ass_dim);
      int a   = a_i + a_j * d_number_y;

      // What pin within the assembly?
      int p_i = i % ass_dim;
      int p_j = j % ass_dim;
      int p   = p_i + p_j * ass_dim;

      // Map to core.
      p += a * ass_dim * ass_dim;
      core_pins[p] = pin;
    }
  } // end assembly loop

  // Fine mesh edges
  vec_dbl edges(number_cells_x + 1, 0.0);

  // Temporary maps.
  vec_int core_mat_map(number_cells, 0);
  vec_int core_reg_map(number_cells, 0);
  vec_int core_pin_map(number_cells, 0);
  vec_int core_ass_map(number_cells, 0);

  // Number of pin cells per assembly
  int number_pins_assembly = d_assemblies[0]->dimension();
  number_pins_assembly *= number_pins_assembly;

  int assembly_count = 0;
  int j_save = 0;
  int i_save = 0;

  // Assembly y-loop
  for (int j = 0; j < number_row; j++)
  {
    // Fine mesh y range for this jth coarse mesh.
    int j1 = j_save;
    int j2 = j_save + d_assemblies[0]->mesh()->number_cells_y();

    // Do edges once, since they are shared.
    for (int jj = j1; jj < j2; jj++)
    {
      // All pins have same mesh in a direction.
      edges[jj + 1] = edges[jj] +
        d_assemblies[assembly_map[j*number_row]]->mesh()->dy(jj - j_save);
    }

    // Assembly x-loop
    for (int i = 0; i < number_row; i++)
    {

      // Pin index
      int pin = i + j*number_row;

      // given assembly i and j, assembly pin id, what's core pin id?

      // Fine mesh x range.
      int i1 = i_save;
      int i2 = i_save + d_assemblies[0]->mesh()->number_cells_x();

      // Get the material and region maps for this pin.
      vec_int ass_mat =
        d_assemblies[assembly_map[pin]]->mesh()->mesh_map("MATERIAL");
      vec_int ass_reg =
        d_assemblies[assembly_map[pin]]->mesh()->mesh_map("REGION");
      vec_int ass_pin =
        d_assemblies[assembly_map[pin]]->mesh()->mesh_map("PINS");

      // Assign the values.
      int count = 0;
      for (int jj = j1; jj < j2; jj++)
      {
        for (int ii = i1; ii < i2; ii++)
        {
          int cell = ii + jj*number_cells_x;
          core_mat_map[cell] = ass_mat[count];
          core_reg_map[cell] = ass_reg[count];
          core_pin_map[cell] = ass_pin[count] +
                               assembly_count * number_pins_assembly;
//          core_pin_map[cell] = core_pins[ass_pin[count] +
//                                 (assembly_count) * number_pins_assembly];
          core_ass_map[cell] = assembly_count;
          count++;
        }
      }
      i_save += d_assemblies[0]->mesh()->number_cells_x();
      assembly_count++;
    }
    i_save = 0;
    j_save += d_assemblies[0]->mesh()->number_cells_y();

  }

  // Create my mesh and add mesh maps.
  d_mesh = new Mesh2D(edges, edges, core_mat_map);
  d_mesh->add_mesh_map("REGION",     core_reg_map);
  d_mesh->add_mesh_map("PINS",       core_pin_map);
  d_mesh->add_mesh_map("ASSEMBLIES", core_ass_map);

}

} // end namespace detran_geometry

