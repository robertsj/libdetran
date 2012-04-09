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
#include "Boundary.hh"
#include "MomentToDiscrete.hh"
#include "State.hh"
#include "Sweeper.hh"
#include "SweepSource.hh"

// Utilities
#include "MathUtilities.hh"

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
 * We have implement the standard source iteration method, Livolant
 * acceleration (an extrapolation technique), and solvers that use
 * MATLAB's own GMRES (and other Krylov solvers).
 *
 * Input parameters specific to InnerIteration and derived classes:
 * - inner_max_iters (default: 100)
 * - inner_tolerance (default: 1e-5)
 *
 * \sa SourceIteration, Livolant, GMRESIteration
 */
// ==============================================================================

template <class D>
class InnerIteration
{

public:

  typedef SP<InnerIteration>                    SP_inner;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename Boundary<D>::SP_boundary     SP_boundary;
  typedef typename MomentToDiscrete<D>::SP_MtoD SP_MtoD;
  // source typedefs
  typedef ExternalSource::SP_source 		    SP_externalsource;
  typedef FissionSource::SP_source 			    SP_fissionsource;
  //
  typedef typename Sweeper<D>::SP_sweeper       SP_sweeper;
  typedef typename
      SweepSource<D>::SP_sweepsource            SP_sweepsource;
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
                 SP_material      	  material,
                 SP_quadrature    	  quadrature,
                 SP_boundary    	  boundary,
                 SP_externalsource    q_e,
                 SP_fissionsource 	  q_f);

  /*!
   *  \brief Solve the within group equation.
   */
  virtual void solve(int g) = 0;

  virtual bool is_valid() const
  { /* ... */ }

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
  /// Group we are solving.
  int d_g;
  /// Numer of sweeps
  int d_number_sweeps;
};

template <class D>
InnerIteration<D>::InnerIteration(SP_input          input,
                                  SP_state          state,
                                  SP_mesh           mesh,
                                  SP_material       material,
                                  SP_quadrature     quadrature,
                                  SP_boundary       boundary,
                                  SP_externalsource q_e,
                                  SP_fissionsource  q_f)
  :  d_input(input)
  ,  d_state(state)
  ,  d_mesh(mesh)
  ,  d_material(material)
  ,  d_quadrature(quadrature)
  ,  d_boundary(boundary)
  ,  d_max_iters(100)
  ,  d_tolerance(1e-5)
{
  Require(input);
  Require(state);
  Require(mesh);
  Require(material);
  Require(quadrature);
  Require(boundary);

  // Get relevant input parameters.
  if (input->check("inner_max_iters"))
  {
    d_max_iters = input->get<int>("inner_max_iters");
  }
  else
  {
    d_max_iters = 5;
  }
  if (input->check("inner_tolerance"))
  {
    d_tolerance = input->get<double>("inner_tolerance");
  }
  else
  {
    d_tolerance = 1e-5;
  }

  // Moments-to-Discrete
  SP_MtoD MtoD;
  MtoD = new MomentToDiscrete<D>(1); // 1 moment (0th order)
  MtoD->build(quadrature);

  // Build the sweep source.
  d_sweepsource =
      new SweepSource<D>(state, mesh, quadrature, material, MtoD);

  if (q_f)
  {
    d_sweepsource->set_fission_source(q_f);
  }

  if (q_e)
  {
    // \todo Come up with a nice way to add/remove sources.  Perhaps
    //       also allow names?
    d_sweepsource->set_moment_source(q_e);
  }

  // Sweeper
  d_sweeper = new Sweeper<D>(input, mesh, material,
		                     quadrature, state, boundary,
		                     d_sweepsource);

}

} // end namespace detran

#endif /* INNERITERATION_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIteration.hh
//---------------------------------------------------------------------------//
