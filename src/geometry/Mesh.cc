/*
 * Mesh.cc
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#include <numeric>

#include "Mesh.hh"

namespace detran
{

Mesh::Mesh(int dim,
           vec_int xfm,  vec_int yfm,  vec_int zfm,
           vec_dbl xcme, vec_dbl ycme, vec_dbl zcme)
  : d_dimension(dim)
  , d_xfm(xfm)
  , d_yfm(yfm)
  , d_zfm(zfm)
  , d_xcme(xcme)
  , d_ycme(ycme)
  , d_zcme(zcme)
{
  Require(d_xfm.size() > 0);
  Require(d_yfm.size() > 0);
  Require(d_zfm.size() > 0);
  Require(d_xcme.size() == d_xfm.size()+1);
  Require(d_ycme.size() == d_yfm.size()+1);
  Require(d_zcme.size() == d_zfm.size()+1);

  // Compute numbers of cells.
  d_number_cells_x = std::accumulate(d_xfm.begin(), d_xfm.end(), 0);
  d_number_cells_y = std::accumulate(d_yfm.begin(), d_yfm.end(), 0);
  d_number_cells_z = std::accumulate(d_zfm.begin(), d_zfm.end(), 0);
  d_number_cells   =  d_number_cells_x * d_number_cells_y * d_number_cells_z;

  // Cell widths.
  d_dx.resize(d_number_cells_x, 0.0);
  d_dy.resize(d_number_cells_y, 0.0);
  d_dz.resize(d_number_cells_z, 0.0);

};


/*!
 *
 */
void Mesh::add_mesh_map(std::string map_key, vec_int &mesh_map)
{
  Require(mesh_map.size() > 0);

  // Erase the value associated with the key if it exists.
  mesh_map_type::iterator it;
  it = d_mesh_map.find(map_key);
  if (it != d_mesh_map.end())
    d_mesh_map.erase(it);

  // Add the new value.
  d_mesh_map[map_key] = mesh_map;
}

}
