/*
 * Mesh2D.cc
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#include <numeric>

#include "Mesh2D.hh"

namespace detran
{

Mesh2D::Mesh2D(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int mat_map)
  : Mesh(2, xfm, yfm, vec_int(1, 1), xcme, ycme, vec_dbl(2, 0.0))
{

  // Does the map cover each coarse mesh?
  Assert(mat_map.size() == d_xfm.size()*d_yfm.size());

  // Temporary fine mesh map.
  vec_int tmp_mat_map(d_number_cells, 0);

  // Discretize
  int ph1 = 0; // place holders
  int ph2 = 0;
  for (int i = 0; i < d_xfm.size(); i++)
  {
    int i1 = ph1;
    int i2 = ph1 + d_xfm[i];
    for (int ii = i1; ii < i2; ii++)
      d_dx[ii] = (d_xcme[i+1] - d_xcme[i]) / d_xfm[i];
    for (int j = 0; j < d_yfm.size(); j++)
    {
      int j1 = ph2;
      int j2 = ph2 + d_yfm[j];
      for (int jj = j1; jj < j2; jj++)
      {
        d_dy[jj] = (d_ycme[j+1] - d_ycme[j]) / d_yfm[j];
        for (int ii = i1; ii < i2; ii++)
          tmp_mat_map[i + j*d_xfm.size()];
      }
      ph2 = std::accumulate(d_yfm.begin(), d_yfm.begin()+j, 0);
    }
    ph2 = 0;
    ph1 = std::accumulate(d_xfm.begin(), d_xfm.begin()+i, 0);
  }
  d_mesh_map["MATERIAL"] = tmp_mat_map;
}

//Mesh2D::Mesh2D(vec_dbl xfme, vec_dbl yfme, vec_int mat_map)
//  : Mesh(2,
//         vec_int(xfme.size()-1, 1),
//         vec_int(yfme.size()-1, 1),
//         vec_int(1, 1),
//         xfme,
//         yfme,
//         vec_dbl(2, 0.0))
//{
//  Require(d_xfm.size() > 0);
//  Require(d_yfm.size() > 0);
//  Require(d_xcme.size() == d_xfm.size()+1);
//  Require(d_ycme.size() == d_yfm.size()+1);
//
//  // Compute numbers of cells.
//  d_number_cells_x = std::accumulate(d_xfm.begin(), d_xfm.end(), 0);
//  d_number_cells_y = std::accumulate(d_yfm.begin(), d_yfm.end(), 0);
//  d_number_cells_z = 1;
//  d_number_cells   =  d_number_cells_x * d_number_cells_y;
//
//  // Cell widths.
//  d_dx.resize(d_number_cells_x, 0.0);
//  d_dy.resize(d_number_cells_y, 0.0);
//  d_dz.resize(1,                1.0);
//
//  // Does the map cover each fine mesh?
//  Assert(mat_map.size() == d_number_cells);
//
//  // Discretize
//  for (int i = 0; i < d_number_cells_x; i++)
//    d_dx[i] = (d_xcme[i+1] - d_xcme[i]) / d_xfm[i];
//  for (int j = 0; j < d_number_cells_y; j++)
//    d_dy[j] = (d_ycme[j+1] - d_ycme[j]) / d_yfm[j];
//
//  // Set material map.
//  d_mesh_map["MATERIAL"] = mat_map;
//
//}

void Mesh2D::add_coarse_mesh_map(std::string map_key, vec_int mesh_map)
{
  return;
}


} // end namespace detran
