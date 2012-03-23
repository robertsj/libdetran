//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Mesh.hh
 * \author Jeremy Roberts
 * \brief  Mesh class definition.
 */
//---------------------------------------------------------------------------//


#ifndef MESH_HH_
#define MESH_HH_

#include "Definitions.hh"
#include "DBC.hh"
#include "SP.hh"
#include "Warning.hh"

#include "map"

namespace detran
{

//using detran_utils::vec_dbl;
//using detran_utils::vec_int;

//===========================================================================//
/*!
 * \class Mesh
 * \brief Abstract Cartesian mesh class.
 *
 */
//===========================================================================//
class Mesh : public detran_utils::Object
{

public:

  enum sides
  {
    LEFT,
    RIGHT,
    BOTTOM,
    TOP,
    SOUTH,
    NORTH
  };

  typedef detran_utils::SP<Mesh>            SP_mesh;
  typedef detran_utils::vec_int             vec_int;
  typedef detran_utils::vec_dbl             vec_dbl;
  typedef std::map<std::string, vec_int>    mesh_map_type;

  /*!
   *  \brief Constructor.
   *
   *  \param    dim         Spatial dimension
   *  \param    xfm         Fine meshes per coarse mesh in x dimension.
   *  \param    yfm         Fine meshes per coarse mesh in y dimension.
   *  \param    zfm         Fine meshes per coarse mesh in z dimension.
   *  \param    xcme        Coarse mesh edges x dimension.
   *  \param    ycme        Coarse mesh edges y dimension.
   *  \param    zcme        Coarse mesh edges z dimension.
   */
  Mesh(int dim,
       vec_int xfm,  vec_int yfm,  vec_int zfm,
       vec_dbl xcme, vec_dbl ycme, vec_dbl zcme);

protected:

  Mesh(int dim) : d_dimension(dim) {}

public:

  //------------------------------------------------------------------------//
  // Setters
  //------------------------------------------------------------------------//

  /*!
   * \brief  Add map of coarse mesh integer properties.
   *
   * This is an easy way to set mesh properties for meshes based on
   * simple coarse mesh regions.  Since the mapping of the coarse mesh
   * to fine mesh is dimension-dependent, we leave it as pure virtual.
   *
   * \param  map_key   String description of map.
   * \param  mesh_map  Logically multi-dimensional map as 1-d vector.
   */
  virtual void add_coarse_mesh_map(std::string map_key, vec_int mesh_map) = 0;

  /*!
   * \brief  Add map of fine mesh integer properties.
   *
   * This adds properties for fine meshes directly, and so is meanty for
   * use with higher level mesh construction, e.g. pin cells, where
   * assignment is not possible by simple coarse mesh bounds.
   *
   * \note If the key exists, this function overwrites the map.
   *
   * \param  map_key   String description of map.
   * \param  mesh_map  Logically multi-dimensional map as 1-d vector.
   */
  void add_mesh_map(std::string map_key, vec_int mesh_map);

  //------------------------------------------------------------------------//
  // Getters
  //------------------------------------------------------------------------//

  int number_cells()
  {
    return d_number_cells;
  }

  int number_cells_x()
  {
    return d_number_cells_x;
  }

  int number_cells_y()
  {
    return d_number_cells_y;
  }

  int number_cells_z()
  {
    return d_number_cells_z;
  }

  double dx(int i)
  {
    Require (i >= 0);
    Require (i < d_number_cells_x);
    return d_dx[i];
  }

  double dy(int j)
  {
    Require (j >= 0);
    Require (j < d_number_cells_y);
    return d_dy[j];
  }

  double dz(int k)
  {
    Require (k >= 0);
    Require (k < d_number_cells_z);
    return d_dz[k];
  }

  int dimension()
  {
    return d_dimension;
  }

  /*!
   *
   * \brief   Returns the cardinal index for i, j, and k
   *
   * \param   i  Index along x axis.
   * \param   j  Index along y axis.
   * \param   k  Index along z axis.
   * \return     Cardinal index.
   */
  virtual int index(int i, int j = 0, int k = 0) = 0;

  /*!
   * \brief  Get map of fine mesh integer properties.
   *
   * This adds properties for fine meshes directly, and so is meanty for
   * use with higher level mesh construction, e.g. pin cells, where
   * assignment is not possible by simple coarse mesh bounds.
   *
   * \param   m  Logically multi-dimensional map as 1-d vector.
   */
  const vec_int& mesh_map(std::string map_key);

  /// Unimplemented DBC function.
  virtual bool is_valid() const {};

protected:

  /// x fine meshes in each x coarse mesh
  vec_int d_xfm;

  /// y fine meshes in each y coarse mesh
  vec_int d_yfm;

  /// z fine meshes in each y coarse mesh
  vec_int d_zfm;

  /// x coarse mesh edges
  vec_dbl d_xcme;

  /// y coarse mesh edges
  vec_dbl d_ycme;

  /// z coarse mesh edges
  vec_dbl d_zcme;

  /// x widths
  vec_dbl d_dx;

  /// y widths
  vec_dbl d_dy;

  /// z widths
  vec_dbl d_dz;

  /// Total number of cells
  int d_number_cells;

  /// Number of cells in x direction
  int d_number_cells_x;

  /// Number of cells in y direction
  int d_number_cells_y;

  /// Number of cells in y direction
  int d_number_cells_z;

  /*!
   *  Map container containing a key describing a mesh property and a fine
   *  mesh map defining the property in each cell.  These properties
   *  include materials, coarse mesh regions (pins, assembly, fuel,
   *  moderator, etc.), and anything else the user wants to edit.
   */
  mesh_map_type d_mesh_map;

  /// Flag indicating I'm meshed.
  bool d_meshed;

  int d_dimension;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//#include "Mesh.i.hh"

#endif /* MESH_HH_ */
