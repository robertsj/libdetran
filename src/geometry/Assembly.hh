//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Assembly.hh
 *  @brief Assembly class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_ASSEMBLY_HH_
#define detran_geometry_ASSEMBLY_HH_

#include "Mesh2D.hh"
#include "PinCell.hh"

namespace detran_geometry
{

/**
 *  @class Assembly
 *  @brief Simple square assembly.
 *
 *  Assembly represents a square array of pin cells that
 *  are assumed to have the same meshing.  The meshing is
 *  also assumed to be isotropic (same in x and y) but not
 *  necessarily uniform.
 */
class GEOMETRY_EXPORT Assembly
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Assembly>    SP_assembly;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef PinCell::SP_pincell               SP_pincell;
  typedef std::vector<SP_pincell>           vec_pincell;
  typedef Mesh::vec_dbl                     vec_dbl;
  typedef Mesh::vec_int                     vec_int;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *
   *  @param    pitch       Pin cell pitch (assumed square)
   *  @param    radii       Vector of fuel pin radii (can be zero length)
   *  @param    mat_map     Region material map (cell-center outward)
   *  @param    meshes      Number of evenly-spaced meshes per direction
   */
  Assembly(int dimension, vec_pincell pincells, vec_int pincell_map);

  /**
   *  @brief Constructor.
   *
   *  @param    dimension   Number of pins per row (e.g 17 for 17x17)
   */
  explicit Assembly(int dimension);

  /// SP Constructor
  static SP_assembly Create(int dim)
  {
    SP_assembly p(new Assembly(dim));
    return p;
  }

  /// Return underlying meshed object.
  Mesh::SP_mesh mesh()
  {
    return d_mesh;
  }

  /// Add a pincell
  void add_pincell(SP_pincell pin);

  /// Mesh the assembly.
  void finalize(vec_int pincell_map);

  /// Get dimension
  int dimension() const
  {
    return d_dimension;
  }

  int number_pincells()
  {
    return d_number_pincells;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Meshed object
  Mesh2D::SP_mesh d_mesh;
  /// Dimension, e.g. 17 in 17x17.
  int d_dimension;
  /// Assembly pitch
  //double d_pitch;
  /// Vector of SP pointers to pin cells in the assembly
  vec_pincell d_pincells;
  /// Logically 2-D map of pin cell locations
  vec_int d_pincell_map;
  /// Number of pincells in the assembly
  int d_number_pincells;

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<Assembly>)
GEOMETRY_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<Assembly> >)

} // end namespace detran_geometry

#endif /* detran_geometry_ASSEMBLY_HH_ */
