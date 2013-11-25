//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Mesh.hh
 *  @brief Mesh class definition.
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_MESH_HH_
#define detran_geometry_MESH_HH_

#include "geometry/geometry_export.hh"
#include "geometry/TrackDB.hh"
#include "geometry/Point.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Warning.hh"
#include <cmath>
#include <map>
#ifdef DETRAN_ENABLE_BOOST
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#endif

namespace detran_geometry
{

//---------------------------------------------------------------------------//
/**
 *  @class Mesh
 *  @brief Abstract Cartesian mesh class.
 *
 *  Note, the constructors are protected to forbid direct instantiation of
 *  the Mesh class.  Rather, use the dimension-specific subclasses.  We
 *  could use a pure virtual destructor as an alternative.
 */
//---------------------------------------------------------------------------//
class GEOMETRY_EXPORT Mesh
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum SIDES
  {
    WEST, EAST, SOUTH, NORTH, BOTTOM, TOP, END_SIDES
  };

  enum FACE2D
  {
    VERT, HORZ
  };

  enum FACE3D
  {
    YZ, XZ, XY
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Mesh>        SP_mesh;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef std::map<std::string, vec_int>    mesh_map_type;
  typedef detran_utilities::size_t          size_t;
  typedef TrackDB::SP_trackdb               SP_trackdb;

  // Note, these constructors are protected to hide them from the
  // user.  These are to be called by inherited classes.  I keep
  // the constructors at the top for convention.
protected:

  /**
   *  @brief Constructor.
   *
   *  @param    dim         Spatial dimension
   *  @param    xfm         Fine meshes per coarse mesh in x dimension.
   *  @param    yfm         Fine meshes per coarse mesh in y dimension.
   *  @param    zfm         Fine meshes per coarse mesh in z dimension.
   *  @param    xcme        Coarse mesh edges x dimension.
   *  @param    ycme        Coarse mesh edges y dimension.
   *  @param    zcme        Coarse mesh edges z dimension.
   */
  Mesh(size_t dim,
       vec_int xfm,  vec_int yfm,  vec_int zfm,
       vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
       vec_int mat_map);

  /**
   *  @brief Constructor.
   *
   *  @param    dim         Spatial dimension
   *  @param    xfme        Fine mesh edges x dimension.
   *  @param    yfme        Fine mesh edges y dimension.
   *  @param    zfme        Fine mesh edges z dimension.
   */
  Mesh(size_t dim,
       vec_dbl xfme, vec_dbl yfme, vec_dbl zfme,
       vec_int mat_map);

  Mesh(size_t dim) : d_dimension(dim) {}

public:

  /// Virtual destructor
  virtual ~Mesh(){}


  //------------------------------------------------------------------------//
  // Setters
  //------------------------------------------------------------------------//

  /**
   * @brief  Add map of coarse mesh integer properties.
   *
   * This is an easy way to set mesh properties for meshes based on
   * the coarse mesh regions used to create the mesh.  If the key exists,
   * this function overwrites the map.
   *
   * @param  map_key   String description of map.
   * @param  mesh_map  Logically multi-dimensional map as 1-d vector.
   */
  void add_coarse_mesh_map(std::string map_key, vec_int mesh_map);

  /**
   * @brief  Add map of fine mesh integer properties.
   *
   * If the key exists, this function overwrites the map.
   *
   * @param  map_key   String description of map.
   * @param  mesh_map  Logically multi-dimensional map as 1-d vector.
   */
  void add_mesh_map(std::string map_key, vec_int mesh_map);

  //------------------------------------------------------------------------//
  // Getters
  //------------------------------------------------------------------------//

  /// Return total number of cells.
  size_t number_cells() const;
  /// Return number of cells in specified dimension.
  size_t number_cells(size_t dim) const;
  /// Return number of cells along x axis
  size_t number_cells_x() const;
  /// Return number of cells along y axis
  size_t number_cells_y() const;
  /// Return number of cells along z axis
  size_t number_cells_z() const;

  /// Get cell width in a specified dimension
  double width(size_t dim, size_t ijk) const;
  /// Get cell width along x axis
  double dx(size_t i) const;
  /// Get cell width along y axis
  double dy(size_t j) const;
  /// Get cell width along z axis
  double dz(size_t k) const;

  /// Get vector of widths along x axis
  const vec_dbl& dx() const;
  /// Get vector of widths along y axis
  const vec_dbl& dy() const;
  /// Get vector of widths along z axis
  const vec_dbl& dz() const;

  /// Get the cell volume
  double volume(size_t cell) const;

  /// Get domain width along x axis
  double total_width_x() const;
  /// Get domain width along y axis
  double total_width_y() const;
  /// Get domain width along z axis
  double total_width_z() const;

  /// Get the mesh dimension
  size_t dimension() const;

  /**
   * @brief   Returns the cardinal index for i, j, and k
   * @param   i  Index along x axis.
   * @param   j  Index along y axis.
   * @param   k  Index along z axis.
   * @return     Cardinal index.
   */
  size_t index(size_t i, size_t j = 0, size_t k = 0);

  /**
   *  @brief   Returns the x index given cardinal index
   *  @param   cell  Cardinal index.
   *  @return        Index along x axis.
   */
  size_t cell_to_i(size_t cell) const;

  /**
   *  @brief   Returns the y index given cardinal index
   *  @param   cell  Cardinal index.
   *  @return        Index along y axis.
   */
  size_t cell_to_j(size_t cell) const;

  /**
   *  @brief   Returns the z index given cardinal index
   *  @param   cell  Cardinal index.
   *  @return        Index along z axis.
   */
  size_t cell_to_k(size_t cell) const;

  /// Find the cell containing a point
  int find_cell(Point p);

  /// Check if fine mesh map exists.
  bool mesh_map_exists(std::string map_key);

  /**
   *  @brief  Get map of fine mesh integer properties.
   *
   *  This adds properties for fine meshes directly, and so is meant for
   *  use with higher level mesh construction, e.g. pin cells, where
   *  assignment is not possible by simple coarse mesh bounds.
   *
   *  @param   m  Logically multi-dimensional map as 1-d vector.
   */
  const vec_int& mesh_map(std::string map_key);

  /// Return a const reference to the full map (useful for IO)
  const mesh_map_type& get_mesh_map() const;

  /// Display some key features
  void display() const;

  void set_tracks(SP_trackdb tracks)
  {
    Require(tracks);
    d_tracks = tracks;
  }

protected:

  /// Mesh spatial dimension
  size_t d_dimension;
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
  /// total x width
  double d_total_width_x;
  /// total y width
  double d_total_width_y;
  /// total z width
  double d_total_width_z;
  /// Total number of cells
  size_t d_number_cells;
  /// Number of cells in x direction
  size_t d_number_cells_x;
  /// Number of cells in y direction
  size_t d_number_cells_y;
  /// Number of cells in y direction
  size_t d_number_cells_z;
  /// Track database
  SP_trackdb d_tracks;

  /**
   *  Map container containing a key describing a mesh property and a fine
   *  mesh map defining the property in each cell.  These properties
   *  include materials, coarse mesh regions (pins, assembly, fuel,
   *  moderator, etc.), and anything else the user wants to edit.
   */
  mesh_map_type d_mesh_map;

private:

  /// Common construction tasks.
  void setup();

#ifdef DETRAN_ENABLE_BOOST

  /// Default constructor.  This is required for serialization and
  /// is of no use otherwise.
  Mesh(){}

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_dimension;
    ar & d_xfm;
    ar & d_yfm;
    ar & d_zfm;
    ar & d_xcme;
    ar & d_ycme;
    ar & d_zcme;
    ar & d_dx;
    ar & d_dy;
    ar & d_dz;
    ar & d_total_width_x;
    ar & d_total_width_y;
    ar & d_total_width_z;
    ar & d_number_cells;
    ar & d_number_cells_x;
    ar & d_number_cells_y;
    ar & d_number_cells_z;
    ar & d_mesh_map;
  }

#endif

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<Mesh>)

} // end namespace detran_geometry

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Mesh.i.hh"

#endif /* detran_geometry_MESH_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh.hh
//---------------------------------------------------------------------------//
