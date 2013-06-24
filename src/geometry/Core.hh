//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   Core.hh
 *  @brief  Core class definition
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_CORE_HH_
#define detran_geometry_CORE_HH_

#include "Mesh2D.hh"
#include "Assembly.hh"

namespace detran_geometry
{

/**
 *  @class Core
 *  @brief Simple 2-D core of assemblies.
 */
class GEOMETRY_EXPORT Core
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Core>    SP_core;
  typedef Mesh::SP_mesh                 SP_mesh;
  typedef Assembly::SP_assembly         SP_assembly;
  typedef std::vector<SP_assembly>      vec_assembly;
  typedef Mesh::vec_int                 vec_int;
  typedef Mesh::vec_dbl                 vec_dbl;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *
   *  @param    pitch       Pin cell pitch (assumed square)
   *  @param    radii       Vector of fuel pin radii (can be zero length)
   *  @param    mat_map     Region material map (cell-center outward)
   *  @param    meshes      Number of evenly-spaced meshes per direction
   */
  Core(int dimension, vec_assembly assemblies, vec_int assembly_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    dimension   Number of pins per row (e.g 17 for 17x17)
   */
  explicit Core(int dimension);

  /// SP Constructor
  static SP_core
  Create(int dimension)
  {
    SP_core p(new Core(dimension));
    return p;
  }

  /// Return underlying meshed object.
  Mesh::SP_mesh mesh()
  {
    return d_mesh;
  }

  /// Add an assembly.
  void add_assembly(SP_assembly assembly);

  vec_assembly assemblies() const
  {
    return d_assemblies;
  }

  /// Get dimension
  int dimension(const size_t dim = 0) const
  {
    if (dim == 0) return d_number_x;
    return d_number_y;
  }

  /// Mesh the assembly.
  void finalize(vec_int assembly_map);

  /// Pincell index.
  int pincell_index(int i, int j)
  {
    return i + j * d_number_x;
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Meshed object
  Mesh2D::SP_mesh d_mesh;
  /// Dimension, e.g. 17 in 17x17. Number assemblies along one dimension.
  int d_number_x;
  int d_number_y;
  /// Vector of SP pointers to pin cells in the assembly
  vec_assembly d_assemblies;
  /// Logically 2-D map of pin cell locations
  vec_int d_assembly_map;
  /// Number of assemblies in the core
  int d_number_assemblies;
  /// Number of pincells in the core
  int d_number_pincells;

};

} // end namespace detran_geometry

#endif /* detran_geometry_CORE_HH_ */
