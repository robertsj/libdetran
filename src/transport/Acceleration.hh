//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Acceleration.hh
 *  @brief  Acceleration class definition
 *  @author Jeremy Roberts
 *  @date   Aug 8, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_ACCELERATION_HH_
#define detran_ACCELERATION_HH_

#include "State.hh"
#include "angle/Quadrature.hh"
#include "discretization/Equation.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

namespace detran
{

/**
 *  @class Acceleration
 *  @brief Base class for coarse mesh acceleration schemes
 *
 *  All anticipated acceleration schemes have several shared
 *  features.  These include
 *    + Computing coarse mesh reaction rates
 *    + Computing functions of coarse mesh boundary fluxes
 *    + Solving a low order equation on the coarse mesh
 */
template <class D>
class Acceleration
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Acceleration>          SP_acceleration;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef detran_angle::Quadrature::SP_quadrature     SP_quadrature;
  typedef State::SP_state                             SP_state;
  typedef typename EquationTraits<D>::face_flux_type  face_flux_type;

  /**
   *  \brief Constructor
   *
   *  @param mesh       Mesh smart pointer
   *  @param material   Material smart pointer
   *  @param quadrature Quadrature smart pointer
   */
  Acceleration(SP_mesh mesh, SP_material material, SP_quadrature quadrature);

  /// Virtual destructor
  virtual ~Acceleration(){}

  /**
   *  @brief Create acceleration mesh given coarseness level and other setup.
   *
   *  By default, this initializes the coarse mesh by
   *  assigning a desired number
   *  of fine meshes per coarse mesh.  Extra meshes are
   *  assigned by round-robin addition until all are assigned.
   *
   *  Clients may re-implement this to do more than just coarsen (e.g.
   *  allocations).
   *
   *  @param level  Desired number of fine meshes per coarse mesh
   */
  virtual void initialize(int level) = 0;

  /*!
   *  \brief Add contribution to an arbitrary function of the coarse
   *         mesh edge flux.
   *
   *  @param  i   x mesh index
   *  @param  j   y mesh index
   *  @param  k   z mesh index
   *  @param  o   octant
   *  @param  a   angle within octant
   *  @param  psi edge angular flux
   */
  virtual void tally(int i, int j, int k, int o, int a, face_flux_type psi) = 0;

  /// Reset for a new sweep.
  virtual void reset() = 0;

  void set_group(int g)
  {
    b_g = g;
  }

  /*!
   *  \brief Homogenize the material data
   *
   *  This function takes the current state vector and homogenizes the
   *  group constants via flux-weighting.
   *
   *  @param state  The current state vector
   */
  //void homogenize(SP_state state, int group);

  /*!
   *  \brief Get the coarse mesh index for a fine mesh
   *  @param  ijk fine mesh index
   *  @param  dim dimension of index
   *  \return     coarse mesh index
   */
  int fine_to_coarse(int ijk, int dim) const;

  /// Return the actual mesh
  SP_mesh get_mesh() const
  {
    return b_mesh;
  }

  /// Return the coarse mesh
  SP_mesh get_coarse_mesh() const
  {
    return b_coarse_mesh;
  }

  /// Return the actual material
  SP_material get_material()
  {
    return b_material;
  }

  /// Return the quadrature
  SP_material get_quadrature()
  {
    return b_quadrature;
  }

  /// Check if the object is in a valid state.
  bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  /// \{

  /// The fine mesh
  SP_mesh b_mesh;

  /// The coarse mesh
  SP_mesh b_coarse_mesh;

  /// Fine mesh material
  SP_material b_material;

  /// Quadrature
  SP_quadrature b_quadrature;

  /// Coarseness level
  int b_level;

  /// Group
  int b_g;

  /// \}

  /// \name Implementation
  /// \{

  /*!
   *  \brief Create the coarse mesh for a given level.
   *  @param  level   Desired number of fine cells per coarse cell
   */
  void coarsen(int level);

  /*!
   *  \brief Check the outgoing edge of a fine mesh cell is on a coarse
   *         mesh boundary.
   *  @param  i   x fine mesh index
   *  @param  j   y fine mesh index
   *  @param  k   z fine mesh index
   *  @param  o   octant index
   */
  bool on_coarse_boundary(int i, int j, int k, int o) const;

  /// \}



};

} // end namespace detran

// Inline member definitions
#include "Acceleration.i.hh"

#endif /* ACCELERATION_HH_ */
