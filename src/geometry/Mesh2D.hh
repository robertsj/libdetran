/*
 * Mesh2D.hh
 *
 *  Created on: Mar 20, 2012
 *      Author: robertsj
 */

#ifndef MESH2D_HH_
#define MESH2D_HH_

#include "map"

#include "Mesh.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class Mesh2D
 * \brief Two-dimensional Cartesian mesh.
 */
//===========================================================================//
class Mesh2D : public Mesh
{

public:

  /*!
   *  \brief Constructor.
   *
   *  \param    xfm         Fine meshes per coarse mesh in x dimension.
   *  \param    yfm         Fine meshes per coarse mesh in y dimension.
   *  \param    xcme        Coarse mesh edges x dimension.
   *  \param    ycme        Coarse mesh edges y dimension.
   *  \param    mat_map     Coarse mesh material map.
   */
  Mesh2D(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int &mat_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    xfme        Fine mesh edges x dimension.
   *  \param    yfme        Fine mesh edges y dimension.
   *  \param    mat_map     Fine mesh material map.
   */
  Mesh2D(vec_dbl xfme, vec_dbl yfme, vec_int &mat_map);

protected:

  /*!
   *   We keep this as an option in the event inherited meshes need
   *   more flexibility.
   */
  Mesh2D() : Mesh(2) {}

public:

  //------------------------------------------------------------------------//
  // Setters
  //------------------------------------------------------------------------//

  /*!
   * \brief  Add map of coarse mesh integer properties.
   *
   * This is an easy way to set mesh properties for meshes based on
   * simple coarse mesh regions.
   *
   * \param  map_key   String description of map.
   * \param  mesh_map  Logically multi-dimensional map as 1-d vector.
   */
  void add_coarse_mesh_map(std::string map_key, vec_int &mesh_map);

  /*!
   *
   * \brief   Returns the cardinal index for i, j, and k
   *
   * \param   i  Index along x axis.
   * \param   j  Index along y axis.
   * \param   k  Index along z axis.
   * \return     Cardinal index.
   */
  int index(int i, int j = 0, int k = 0)
  {
    Require(i >= 0);
    Require(i < d_number_cells_x);
    Require(j >= 0);
    Require(j < d_number_cells_y);
    return i + j * d_number_cells_x;
  }


};

} // end namespace detran

#endif /* MESH2D_HH_ */
