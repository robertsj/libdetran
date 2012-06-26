//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER_HH_
#define SWEEPER_HH_

// Config
#include "detran_config.h"

// Detran
#include "BoundaryBase.hh"
#include "Equation.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"
#include "SweepSource.hh"

// Detran Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "InputDB.hh"
#include "Profiler.hh"
#include "SP.hh"

// System
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

  typedef SP<Sweeper>                       SP_sweeper;
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  //
  typedef BoundaryBase<D>                   Boundary_T;
  //
  typedef typename
      SweepSource<D>::SP_sweepsource        SP_sweepsource;
  //
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

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
  virtual void setup_group(int g) = 0;

  /// Allows the psi update to occur whenever needed
  void set_update_psi(bool v)
  {
    d_update_psi = v;
  }

  /// Switch on-the-fly boundary updates on or off
  void set_update_boundary(bool v)
  {
    d_update_boundary = v;
  }

  bool update_boundary() const
  {
    return d_update_boundary;
  }

  int number_sweeps()
  {
    return d_number_sweeps;
  }

  /// Set an acceleration scheme.
//  void set_acceleration(SP_acceleration acceleration)
//  {
//    d_acceleration = acceleration;
//  }

  bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  /// \{

  /// Input database
  SP_input d_input;

  /// Material
  SP_material d_material;

  /// Mesh
  SP_mesh d_mesh;

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
  int d_g;

  /// Update the angular flux?
  bool d_update_psi;

  /// Match incident/outgoing side with octant
  vec3_int d_face_index;

  /// Adjoint problem?
  bool d_adjoint;

  /// Count the sweeps.
  int d_number_sweeps;

  /// Update the boundary on the fly?  Can't be used for Krylov.
  bool d_update_boundary;

  /// \}

  /// \name Implementation
  /// \{

  /// Allocate template-specific items.
  void setup();

  /// Mesh sweeper indices. \todo Allow adjoint.
  inline int index(int o, int dim, int ijk);

  /// \}

};

} // end namespace detran

#include "Sweeper.t.hh"

#endif /* SWEEPER_HH_ */
