//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Assembly.cc
 *  @brief Assembly class member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Assembly.hh"
#include <cmath>
#include <iostream>

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Assembly::Assembly(const size_t nx, const size_t ny)
  : d_number_x(nx)
  , d_number_y(ny)
{
  Require(d_number_x > 0);
  if (d_number_y == 0) d_number_y = d_number_x;
}

//----------------------------------------------------------------------------//
Assembly::SP_assembly Assembly::Create(const size_t nx, const size_t ny)
{
  SP_assembly p(new Assembly(nx, ny));
  return p;
}

//----------------------------------------------------------------------------//
Assembly::SP_mesh Assembly::mesh()
{
  return d_mesh;
}

//----------------------------------------------------------------------------//
void Assembly::add_pincell(SP_pincell pin)
{
  Require(pin);
  d_pincells.push_back(pin);
}

//----------------------------------------------------------------------------//
void Assembly::set_pincell_map(const vec_int &pincell_map)
{
  // Verify the pin cell map is correct size.  This does not check to see
  // if pins match up.
  Insist(pincell_map.size() == d_number_x * d_number_y,
         "Pincell map is the wrong size.");
  d_pincell_map = pincell_map;


}

Assembly::vec_int Assembly::pincell_map() const
{
  return d_pincell_map;
}

//----------------------------------------------------------------------------//
void Assembly::finalize(vec_int pincell_map)
{
  using std::cout;
  using std::endl;

  // Verify that pin cells are consistent, using first
  // pin as the reference.
  for (size_t i = 1; i < d_pincells.size(); i++)
  {
    Insist(d_pincells[0]->mesh()->number_cells_x() ==
           d_pincells[i]->mesh()->number_cells_x(),
           "X fine mesh inconsistent.");
    Insist(d_pincells[0]->mesh()->number_cells_y() ==
           d_pincells[i]->mesh()->number_cells_y(),
           "Y fine mesh inconsistent.");
  }

  // Verify the pin cell map is correct.
  Insist(pincell_map.size() == d_number_x*d_number_y,
         "Pincell map is the wrong size.");
  d_pincell_map = pincell_map;

  size_t number_pincells = d_number_x * d_number_y;

  // Set number of cells.  This *assumes* all pins have the same meshing,
  // which was already checked above.
  int number_pins_row = std::sqrt(pincell_map.size());
  Assert(pincell_map.size() == number_pins_row*number_pins_row);
  int number_cells_x =
    (d_pincells[0]->mesh())->number_cells_x() * number_pins_row;
  int number_cells_y =
    (d_pincells[0]->mesh())->number_cells_x() * number_pins_row;
  int number_cells = number_cells_x * number_cells_y;

  // Fine mesh edges
  vec_dbl edges(number_cells_x + 1, 0.0);

  // Temporary fine mesh maps.
  vec_int ass_mat_map(number_cells, 0);
  vec_int ass_reg_map(number_cells, 0);
  vec_int ass_pin_map(number_cells, 0);

  int pin_count = 0;
  int j_save = 0;
  int i_save = 0;

  for (int j = 0; j < number_pins_row; j++)
  {
    // Fine mesh y range for this jth coarse mesh.
    int j1 = j_save;
    int j2 = j_save + d_pincells[0]->mesh()->number_cells_y();

    // Do edges once, since they are shared by x and y.
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
      vec_int pin_mat_map =
        d_pincells[pincell_map[pin]]->mesh()->mesh_map("MATERIAL");
      vec_int pin_reg_map =
        d_pincells[pincell_map[pin]]->mesh()->mesh_map("REGION");

      // Is this a fuel pin?  If not, assign -1 to the pin map.  Then, 0..N are
      // the actual fuel pins.
      int pin_id = -1;
      if (1==1)//d_pincells[pincell_map[pin]]->is_fuel())
      {
        pin_id = pin_count;
        pin_count++;
      }

      // Assign the values.
      int count = 0;
      for (int jj = j1; jj < j2; jj++)
      {
        for (int ii = i1; ii < i2; ii++)
        {
          int cell = ii + jj*number_cells_x;
          ass_mat_map[cell] = pin_mat_map[count];
          ass_reg_map[cell] = pin_reg_map[count];
          ass_pin_map[cell] = pin_id;
          count++;
        }
      }

      i_save += d_pincells[0]->mesh()->number_cells_x();
    }

    i_save = 0;
    j_save += d_pincells[0]->mesh()->number_cells_y();
  }

  // Create my mesh.
  d_mesh = new Mesh2D(edges, edges, ass_mat_map);
  // Add maps.
  d_mesh->add_mesh_map("REGION", ass_reg_map);
  // Assigns unique edit region for each pin in an assembly.
  d_mesh->add_mesh_map("PINS", ass_pin_map);

}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Assembly.hh
//----------------------------------------------------------------------------//
