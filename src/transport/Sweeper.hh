//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper class definition.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER_HH_
#define SWEEPER_HH_

#include "detran_config.h"
#include "CurrentTally.hh"
#include "State.hh"
#include "SweepSource.hh"
#include "angle/Quadrature.hh"
#include "boundary/BoundaryBase.hh"
#include "discretization/Equation.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Sweeper
 * \brief Sweeper for discrete ordinates problems.
 *
 * The within-group transport equation is
 * \f[
 *      \mathbf{L}\psi = Q \, ,
 * \f]
 * where \f$ \mathbf{L} \f$ is the streaming and collision operator and
 * \f$ Q \f$ is a discrete representation of all source contributions.
 *
 * To invert the operator \f$ \mathbf{L} \f$, we "sweep" over the mesh for all
 * angles,
 * which gives us updated angular fluxes in each cell.  Actually, the
 * flux *moments* are updated, while the discrete angular flux is
 * optionally stored.
 *
 * Relevant input database entries:
 *   - store_angular_flux [int]
 *   - equation [string]
 *
 */
//---------------------------------------------------------------------------//
template <class D>
class Sweeper
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Sweeper>             SP_sweeper;
  typedef State::SP_state                           SP_state;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh                     Mesh;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef BoundaryBase<D>                           Boundary_T;
  typedef typename SweepSource<D>::SP_sweepsource   SP_sweepsource;
  typedef State::moments_type                       moments_type;
  typedef State::angular_flux_type                  angular_flux_type;
  typedef CurrentTally<D>                           CurrentTally_T;
  typedef typename CurrentTally_T::SP_currenttally  SP_currenttally;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::vec2_int                vec2_int;
  typedef detran_utilities::vec3_int                vec3_int;
  typedef detran_utilities::vec_size_t              vec_size_t;
  typedef detran_utilities::vec2_size_t             vec2_size_t;
  typedef detran_utilities::vec3_size_t             vec3_size_t;
  typedef detran_utilities::size_t                  size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature.
   *  \param    state       State vectors.
   */
  Sweeper(SP_input input,
          SP_mesh mesh,
          SP_material material,
          SP_quadrature quadrature,
          SP_state state,
          SP_sweepsource sweepsource);

  /// Virtual destructor
  virtual ~Sweeper(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SWEEPERS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Sweep over all angles and space.
   *
   *  Note, if the angular flux is to be updated,
   *  it is done directly to via State.  Having
   *  sweep take the flux as an explicit argument
   *  allows various input types (e.g. Krylov
   *  vectors) without having to go through State.
   *
   *  Note, here the source is limited to a moments-based
   *  one.  We may want a discrete source later; that can
   *  be added with a third argument of psi type.
   *
   */
  virtual void sweep(moments_type &phi) = 0;

  /// Setup the equations for the group
  virtual void setup_group(const size_t g) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /// Allows the psi update to occur whenever needed
  void set_update_psi(const bool v)
  {
    d_update_psi = v;
  }

  /// Switch on-the-fly boundary updates on or off
  void set_update_boundary(const bool v)
  {
    d_update_boundary = v;
  }

  bool update_boundary() const
  {
    return d_update_boundary;
  }

  size_t number_sweeps() const
  {
    return d_number_sweeps;
  }

  /// Set a current tally.
  void set_current(SP_currenttally current)
  {
    Require(current);
    d_current = current;
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Input database
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Material
  SP_material d_material;

  /// Angular quadrature
  SP_quadrature d_quadrature;

  /// State vectors
  SP_state d_state;

  /// Boundary
  //SP_boundary d_boundary;

  /// Sweep source
  SP_sweepsource d_sweepsource;

  /// Acceleration
  //SP_acceleration d_acceleration;

  /// Current group
  size_t d_g;

  /// Update the angular flux?
  bool d_update_psi;

  /// Match incident/outgoing side with octant
  vec3_size_t d_face_index;

  /// Adjoint problem?
  bool d_adjoint;

  /// Count the sweeps.
  size_t d_number_sweeps;

  /// Update the boundary on the fly?  Can't be used for Krylov.
  bool d_update_boundary;

  /// Current tally
  SP_currenttally d_current;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Allocate template-specific items.
  void setup();

  /// Mesh sweeper indices. \todo Allow adjoint.
  inline size_t index(const size_t o, const size_t dim, const size_t ijk);

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "Sweeper.t.hh"

#endif /* SWEEPER_HH_ */
