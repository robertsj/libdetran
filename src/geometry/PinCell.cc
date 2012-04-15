/*
 * PinCell.cc
 *
 *  Created on: Apr 14, 2012
 *      Author: robertsj
 */

// Detran
#include "PinCell.hh"

// System
#include <cmath>

namespace detran
{

PinCell::PinCell(double pitch, vec_dbl radii, vec_int mat_map, int number_meshes)
  : d_pitch(pitch)
  , d_radii(radii)
{
  Require(d_pitch > 0.0);
  Require(number_meshes > 0);
  Require(mat_map.size() == 1 + radii.size());
  //
  double width = d_pitch / number_meshes;
  d_dx.resize(number_meshes, width);
  d_dy.resize(number_meshes, width);
  d_dz.resize(1, 1.0);
  d_number_cells_x = number_meshes;
  d_number_cells_y = number_meshes;
  d_number_cells_z = 1;
  d_number_cells   = number_meshes * number_meshes;

  // Temporary fine mesh material map
  vec_int tmp_mat_map(d_number_cells, 1);
  vec_int tmp_reg_map(d_number_cells, 1);

  //double volume = d_number_regions, 1);
  for (int j = 0; j < number_meshes; j++)
  {
    for (int i = 0; i < number_meshes; i++)
    {
      // Which region am I in?
      int r = find_region(i, j);
      // Which material does this region have?
      int m = mat_map[r];
      // Assign the values.
      int cell = index(i, j, 0);
      tmp_mat_map[cell] = m;
      tmp_reg_map[cell] = r;
    }
  }
  // Add maps
  add_mesh_map("MATERIAL", tmp_mat_map);
  add_mesh_map("REGION", tmp_reg_map);
}


int PinCell::find_region(int i, int j)
{
  double x = (i + 0.5) * d_dx[0];
  double y = (j + 0.5) * d_dx[0];
  // Loop through the radii. If I'm in there, that's where I live.
  int r = d_radii.size(); // start off in outer region.
  double hp = 0.5 * d_pitch;
  for (int p = 0; p < d_radii.size(); p++)
  {
    if (std::sqrt((x - hp) * (x - hp) + (y - hp) * (y - hp)) < d_radii[p])
    {
      r = p;
      break;
    }
  }
  return r;
}

} // end namespace detran



