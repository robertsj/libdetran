/*
 * Core.hh
 *
 *  Created on: Apr 16, 2012
 *      Author: robertsj
 */

#ifndef CORE_HH_
#define CORE_HH_

#include "Mesh2D.hh"
#include "Assembly.hh"

namespace detran
{

/*!
 *  \class Core
 *  \brief Simple 2-D core of assemblies.
 */
class Core : public Object
{

public:

  typedef SP<Core>                  SP_core;
  typedef Mesh::SP_mesh             SP_mesh;
  typedef Assembly::SP_assembly     SP_assembly;
  typedef std::vector<SP_assembly>  vec_assembly;

  /*!
   *  \brief Constructor.
   *
   *  \param    pitch       Pin cell pitch (assumed square)
   *  \param    radii       Vector of fuel pin radii (can be zero length)
   *  \param    mat_map     Region material map (cell-center outward)
   *  \param    meshes      Number of evenly-spaced meshes per direction
   */
  Core(int dimension, vec_assembly assemblies, vec_int assembly_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    dimension   Number of pins per row (e.g 17 for 17x17)
   */
  explicit Core(int dimension);

  /// SP Constructor
  static SP<Core> Create(int dimension)
  {
    SP_core p;
    p = new Core(dimension);
    return p;
  }

  /// Return underlying meshed object.
  Mesh::SP_mesh mesh()
  {
    return d_mesh;
  }

  /// Get const reference to my mesh.
  const Mesh2D& mesh_ref() const
  {
    return *d_mesh;
  }

  /// Add an assembly.
  void add_assembly(SP_assembly assembly);

  /// Mesh the assembly.
  void finalize(vec_int assembly_map);

  bool is_valid() const
  {
    return true;
  }

private:

  /// Meshed object
  Mesh2D::SP_mesh d_mesh;

  /// Dimension, e.g. 17 in 17x17.
  int d_dimension;

  /// Vector of SP pointers to pin cells in the assembly
  vec_assembly d_assemblies;

  /// Logically 2-D map of pin cell locations
  vec_int d_assembly_map;
};

} // end namespace detran

#endif /* CORE_HH_ */
