//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIteration.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerIteration class definition.
 */
//---------------------------------------------------------------------------//

#ifndef INNERITERATION_HH_
#define INNERITERATION_HH_

#include "angle/MomentToDiscrete.hh"
#include "boundary/BoundaryBase.hh"
#include "external_source/ExternalSource.hh"
#include "transport/State.hh"
#include "transport/Sweeper.hh"
#include "transport/Sweeper1D.hh"
#include "transport/Sweeper2D.hh"
#include "transport/Sweeper3D.hh"
#include "transport/Sweeper2DMOC.hh"
#include "transport/SweepSource.hh"
#include "utilities/InputDB.hh"
#include "utilities/MathUtilities.hh"
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
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
//---------------------------------------------------------------------------//

template <class D>
class InnerIteration
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<InnerIteration>        SP_inner;
  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef detran_angle::Quadrature::SP_quadrature     SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_angle::MomentToDiscrete::SP_MtoD     SP_MtoD;
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_externalsource;
  typedef FissionSource::SP_fissionsource 		        SP_fissionsource;
  typedef typename Sweeper<D>::SP_sweeper             SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource     SP_sweepsource;
  typedef State::moments_type                         moments_type;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL WITHIN-GROUP SOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Solve the within group equation.
   */
  virtual void solve(const size_t g) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

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
  size_t d_max_iters;
  /// Convergence tolerance
  double d_tolerance;

  /// Print out flag
  size_t d_print_out;
  /// Interval for print out
  size_t d_print_interval;

  /// Group we are solving.
  size_t d_g;
  /// Number of sweeps
  size_t d_number_sweeps;

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
