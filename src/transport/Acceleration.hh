/*
 * Acceleration.hh
 *
 *  Created on: May 16, 2012
 *      Author: robertsj
 */

#ifndef ACCELERATION_HH_
#define ACCELERATION_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"

// Utilities
#include "DBC.hh"

namespace detran
{

class Acceleration : public Object
{

public:

  /// \name Useful Typedefs
  // \{

  Mesh::SP_mesh SP_mesh;

  Material::SP_material SP_material;

  State::SP_state SP_state;

  // \}

  Acceleration(SP_mesh m, SP_material mat);

  ~Acceleration(){}

  /*!
   *  \brief Create acceleration mesh given coarseness level
   *
   *  Initialize the coarse mesh by assigning a desired number
   *  of fine meshes per coarse mesh.  Extra meshes are
   *  assigned by round-robin addition until all are assigned.
   *
   *  \param level  Desired number of fine meshes per coarse mesh
   */
  void initialize(int level);

  /*!
   *  \brief Homogenize the material data
   *
   */
  void homogenize(SP_state state);

  /*!
   *  \brief Get the coarse mesh index for a fine mesh
   *  \param  ijk fine mesh index
   *  \param  dim dimension of index
   *  \return     coarse mesh index
   */
  int fine_to_coarse(int ijk, int dim);

  /// Return the actual mesh
  SP_mesh mesh()
  {
    return d_mesh;
  }

  /// Return the coarse mesh
  SP_mesh course_mesh()
  {
    return d_coarse_mesh;
  }

  /// Return the actual material
  SP_material material()
  {
    return d_material;
  }

  /// Return the coarse mesh material
  SP_material coarse_material()
  {
    return d_coarse_material;
  }

  bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  //

  /// The fine mesh
  SP_mesh d_mesh;

  /// The coarse mesh
  SP_mesh d_coarse_mesh;

  /// Fine mesh material
  SP_material d_material;

  /// Coarse mesh material
  SP_material d_coarse_material;

  /// Fine-to-coarse maps
  vec_int d_fine_to_coarse_x;
  vec_int d_fine_to_coarse_y;
  vec_int d_fine_to_coarse_z;

};

} // end namespace detran


#endif /* ACCELERATION_HH_ */
