//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIteration.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef INNERITERATION_HH_
#define INNERITERATION_HH_

// Detran
//#include "Acceleration.hh"
#include "BoundaryBase.hh"
#include "MomentToDiscrete.hh"
#include "State.hh"
#include "Sweeper.hh"
#include "Sweeper1D.hh"
#include "Sweeper2D.hh"
#include "Sweeper3D.hh"
#include "Sweeper2DMOC.hh"
#include "SweepSource.hh"

// Utilities
#include "InputDB.hh"
#include "MathUtilities.hh"

// System
#include <string>

namespace detran
{

//===========================================================================//
/*!
 * \class InnerIteration
 * \brief Solve the within-group transport equation.
 *
 * The within-group transport equation in operator form is
 * \f[
 *      \mathbf{L}\psi =
 *      \mathbf{MS}\phi + Q
 * \f]
 * where \f$ \mathbf{L} \f$ is the streaming and collision operator,
 * \mathbf{M} is the moment-to-discrete operator,
 * \mathbf{S} is the scattering operator, and  Q represents any
 * source considered fixed, which includes in-scatter, fission, and
 * external sources.
 *
 * What we are really after is the scalar flux and possibly its higher
 * order moments.  Consequently, we are able to solve a somewhat different
 * problem then the within group transport equation above.  Let us operate
 * on both sides by \f$\mathbf{L}^{-1}\f$ followed
 * by \f$ \mathbf{D}\f$ to get
 * \f[
 *      (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi
 *      = \mathbf{D} \mathbf{L}^{-1} Q \, .
 * \f]
 * Here, \f$\mathbf{D}\f$ is the discrete-to-moment operator, defined
 * such that \f$ \phi = \mathbf{D}\psi \f$.
 *
 * Notice this is nothing but a linear system of the form
 * \f$ \mathbf{A}x = b \f$ where
 * \f[
 *      \mathbf{A} = (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})
 * \f]
 * and
 * \f[
 *      b = \mathbf{D} \mathbf{L}^{-1} Q \, .
 * \f]
 * Moreover, \f$ b\f$ is just the uncollided flux.
 *
 * A nice overview of approaches for this inner iteration is
 * given by Larsen and Morel in <em> Nuclear Computational Science </em>.
 *
 * Input parameters specific to InnerIteration and derived classes:
 * - inner_max_iters        [100]
 * - inner_tolerance        [1e-5]
 * - inner_print_out        [2], 0=never, 1=final, 2=every interval
 * - inner_print_interval   [10]
 *
 * \sa SourceIteration, InnerGMRES
 */
// ==============================================================================

template <class D>
class InnerIteration: public Object
{

public:

  typedef SP<InnerIteration>                    SP_inner;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  typedef typename MomentToDiscrete<D>::SP_MtoD SP_MtoD;
  // source typedefs
  typedef ExternalSource::SP_source 		        SP_externalsource;
  typedef FissionSource::SP_source 			        SP_fissionsource;
  // sweep
  typedef typename Sweeper<D>::SP_sweeper       SP_sweeper;
  typedef typename
      SweepSource<D>::SP_sweepsource            SP_sweepsource;
  // acceleration
//  typedef typename
//      Acceleration<D>::SP_acceleration          SP_acceleration;
  //
  typedef State::moments_type                   moments_type;


  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param external_source   User-defined external source.
   *  \param fission_source    Fission source.
   */
  InnerIteration(SP_input       	  input,
                 SP_state       	  state,
                 SP_mesh        	  mesh,
                 SP_material      	material,
                 SP_quadrature      quadrature,
                 SP_boundary    	  boundary,
                 SP_externalsource  q_e,
                 SP_fissionsource   q_f);


  virtual ~InnerIteration(){}

  /*!
   *  \brief Solve the within group equation.
   */
  virtual void solve(int g) = 0;

  /// \name Setters
  /// \{

  /// Reset the tolerance.
  void set_tolerance(double tol)
  {
    Require(tol > 0.0);
    d_tolerance = tol;
  }

  /// Reset the maximum iterations.
  void set_max_iters(int max_iters)
  {
    Require(max_iters > 0);
    d_max_iters = max_iters;
  }

  /// \}

  /// \name Getters
  /// \{

  SP_mesh get_mesh() const
  {
    return d_mesh;
  }

  SP_state get_state() const
  {
    return d_state;
  }

  SP_sweepsource get_sweepsource() const
  {
    return d_sweepsource;
  }

  SP_sweeper get_sweeper() const
  {
    return d_sweeper;
  }

  int group() const
  {
    return d_g;
  }

  /// \}

  /// Unimplemented DBC function.
  virtual bool is_valid() const
  {
    return true;
  }

protected:

  /// User input.
  SP_input d_input;
  /// State vectors.
  SP_state d_state;
  /// Problem mesh (either Cartesian mesh or MOC tracking)
  SP_mesh d_mesh;
  /// Materials definitions.
  SP_material d_material;
  /// Angular mesh.
  SP_quadrature d_quadrature;
  /// Boundary fluxes.
  SP_boundary d_boundary;

  /// Sweeper over the space-angle domain.
  SP_sweeper d_sweeper;
  /// Sweep source constructor
  SP_sweepsource d_sweepsource;
  /// User-defined external source
  SP_externalsource d_external_source;
  /// Fission source, if used
  SP_fissionsource d_fission_source;

  /// Maximum iterations
  int d_max_iters;
  /// Convergence tolerance
  double d_tolerance;

  /// Print out flag
  int d_print_out;
  /// Interval for print out
  int d_print_interval;

  /// Group we are solving.
  int d_g;
  /// Number of sweeps
  int d_number_sweeps;

  /// Low order acceleration
 // SP_acceleration b_acceleration;

private:

  /// Setup the templated sweeper
  inline bool set_sweep(std::string equation);

};

} // end namespace detran

#include "InnerIteration.t.hh"

#endif /* INNERITERATION_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIteration.hh
//---------------------------------------------------------------------------//
