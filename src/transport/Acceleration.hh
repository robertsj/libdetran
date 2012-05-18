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
#include "Quadrature.hh"
#include "State.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

/*!
 *  \class Acceleration
 *  \brief Base class for coarse mesh acceleration schemes
 *
 *  All anticipated acceleration schemes have several shared
 *  features.  These include
 *    + Requiring reaction rates within a coarse mesh
 *    + Requiring knowledge of the angular flux at coarse
 *      mesh boundaries (to compute net currents, partial
 *      currents, or some other function of the flux)
 *    + Solution of some lower order equation with the
 *      condition that the lower order solution is
 *      equivelent to the homogenized (and converged)
 *      high order solution
 *
 */
class Acceleration : public Object
{

public:

  /// \name Useful Typedefs
  // \{

  typedef SP<Acceleration>          SP_acceleration;
  typedef Mesh::SP_mesh             SP_mesh;
  typedef Material::SP_material     SP_material;
  typedef Quadrature::SP_quadrature SP_quadrature;
  typedef State::SP_state           SP_state;

  // \}

  /*!
   *  \brief Constructor
   *
   *  \param mesh       Mesh smart pointer
   *  \param material   Material smart pointer
   *  \param quadrature Quadrature smart pointer
   */
  Acceleration(SP_mesh mesh, SP_material material, SP_quadrature quadrature);

  /// Virtual destructor
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
  virtual void initialize(int level);

  /*!
   *  \brief Add contribution to an arbitrary function of the coarse
   *         mesh edge flux.
   *
   *  \param  i   x mesh index
   *  \param  j   y mesh index
   *  \param  k   z mesh index
   *  \param  psi edge angular flux
   *  \param  o   octant
   *  \param  a   angle within octant
   *  \param  g   group
   */
  virtual void tally(int i, int j, int k, int o, int a, int g, double psi) = 0;

  /*!
   *  \brief Homogenize the material data
   *
   *  This function takes the current state vector and homogenizes the
   *  group constants via flux-weighting.
   *
   *  \param state  The current state vector
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
  SP_mesh get_mesh()
  {
    return b_mesh;
  }

  /// Return the coarse mesh
  SP_mesh get_coarse_mesh()
  {
    return b_coarse_mesh;
  }

  /// Return the actual material
  SP_material get_material()
  {
    return b_material;
  }

  /// Return the coarse mesh material
  SP_material get_coarse_material()
  {
    return b_coarse_material;
  }

  /// Check if the object is in a valid state.
  bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  //

  /// The fine mesh
  SP_mesh b_mesh;

  /// The coarse mesh
  SP_mesh b_coarse_mesh;

  /// Fine mesh material
  SP_material b_material;

  /// Coarse mesh material
  SP_material b_coarse_material;

  /// Quadrature
  SP_quadrature b_quadrature;

  /// Fine-to-coarse maps
  vec_int b_fine_to_coarse_x;
  vec_int b_fine_to_coarse_y;
  vec_int b_fine_to_coarse_z;

  /// Coarseness level
  int b_level;

};

} // end namespace detran

// Inline member definitions
#include "Acceleration.i.hh"

#endif /* ACCELERATION_HH_ */
