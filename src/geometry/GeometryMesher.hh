//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  GeometryMesher.hh
 *  @brief GeometryMesher class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_GEOMETRYMESHER_HH_
#define detran_geometry_GEOMETRYMESHER_HH_

#include "geometry/Geometry.hh"
#include "geometry/Mesh.hh"

namespace detran_geometry
{

/**
 *  @class GeometryMesher
 *  @brief Methods to discretize a @ref Geometry on a Cartesian mesh
 */
class GEOMETRY_EXPORT GeometryMesher
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATION
  //--------------------------------------------------------------------------//

  enum MESHING_SCHEMES
  {
    CONSTANT, OPTIMIZED, END_MESHING_SCHEMES
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Geometry::SP_geometry             SP_geometry;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_size_t      vec_size_t;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec_size_t      vec_size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Constructor
  GeometryMesher(const double max_spacing,
                 const size_t max_cells);

  /// Meshify the geometry
  void meshify(SP_geometry  geo, const size_t scheme = CONSTANT);

  SP_mesh mesh() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Maximum spacing allowed in any direction
  double d_max_spacing;
  /// Maximum number of cells in any direction
  size_t d_max_cells;
  //@{
  /// Mesh edges
  vec_dbl d_x_edges;
  vec_dbl d_y_edges;
  vec_dbl d_z_edges;
  //@}
  /// Map specifying materials in each region
  vec_size_t d_material_map;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  void meshify_2D(SP_geometry geo, const size_t scheme);

};


} // end namespace detran_geometry

#endif /* detran_geometry_GEOMETRYMESHER_HH_ */
