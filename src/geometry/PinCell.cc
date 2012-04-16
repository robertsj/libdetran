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

// Constructor
PinCell::PinCell(double pitch, vec_dbl radii, vec_int mat_map)
  : d_pitch(pitch)
  , d_radii(radii)
  , d_mat_map(mat_map)
{
  Require(d_pitch > 0.0);
  Require(d_mat_map.size() == 1 + d_radii.size());
}

// Mesh the object
void PinCell::meshify(int number_meshes)
{
  Require(number_meshes > 0);

  // Fine mesh width.
  double width = d_pitch / number_meshes;

  // Fine mesh edges.  Assumes constant spacing.
  vec_dbl edges(number_meshes + 1, 0.0);

  // Number of cells in the meshed pin cell.
  int number_cells = number_meshes*number_meshes;

  // Temporary fine mesh material map
  vec_int tmp_mat_map(number_cells, 1);
  vec_int tmp_reg_map(number_cells, 1);

  //double volume = d_number_regions, 1);
  for (int j = 0; j < number_meshes; j++)
  {
    edges[j+1] = edges[j] + width;
    for (int i = 0; i < number_meshes; i++)
    {
      // Which region am I in?
      int r = find_region(i, j, width);
      // Which material does this region have?
      int m = d_mat_map[r];
      // Assign the values.
      int cell = i + j * number_meshes;
      tmp_mat_map[cell] = m;
      tmp_reg_map[cell] = r;
    }
  }

  // Create my mesh.
  d_mesh = new Mesh2D(edges, edges, tmp_mat_map);

  // Add maps
  d_mesh->add_mesh_map("MATERIAL", tmp_mat_map);
  d_mesh->add_mesh_map("REGION", tmp_reg_map);

}


int PinCell::find_region(int i, int j, double width)
{
  double x = (i + 0.5) * width;
  double y = (j + 0.5) * width;
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



