/*
 * Assembly.hh
 *
 *  Created on: Apr 14, 2012
 *      Author: robertsj
 */

#ifndef ASSEMBLY_HH_
#define ASSEMBLY_HH_

#include "Mesh2D.hh"
#include "PinCell.hh"

namespace detran
{

/*!
 *  \class Assembly
 *  \brief Simple square assembly.
 *
 *  Assembly represents a square array of pin cells that
 *  are assumed to have the same meshing.  The meshing is
 *  also assumed to be isotropic (same in x and y) but not
 *  necessarily uniform.
 *
 */
class Assembly : public Object
{

public:

  typedef SP<Assembly>            SP_assembly;
  typedef Mesh::SP_mesh           SP_mesh;
  typedef PinCell::SP_pincell     SP_pincell;
  typedef std::vector<SP_pincell> vec_pincell;

  /*!
   *  \brief Constructor.
   *
   *  \param    pitch       Pin cell pitch (assumed square)
   *  \param    radii       Vector of fuel pin radii (can be zero length)
   *  \param    mat_map     Region material map (cell-center outward)
   *  \param    meshes      Number of evenly-spaced meshes per direction
   */
  Assembly(int dimension, vec_pincell pincells, vec_int pincell_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    dimension   Number of pins per row (e.g 17 for 17x17)
   */
  explicit Assembly(int dimension);

  /// SP Constructor
  static SP<Assembly> Create(int dim)
  {
    SP_assembly p;
    p = new Assembly(dim);
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

  /// Add a pincell
  void add_pincell(SP_pincell pin);

  /// Mesh the assembly.
  void finalize(vec_int pincell_map);

  /// Get dimension
  int dimension() const
  {
    return d_dimension;
  }

  bool is_valid() const
  {
    return true;
  }

private:

  /// Meshed object
  Mesh2D::SP_mesh d_mesh;

  /// Dimension, e.g. 17 in 17x17.
  int d_dimension;

  /// Assembly pitch
  double d_pitch;

  /// Vector of SP pointers to pin cells in the assembly
  vec_pincell d_pincells;

  /// Logically 2-D map of pin cell locations
  vec_int d_pincell_map;
};

} // end namespace detran

#endif /* ASSEMBLY_HH_ */
