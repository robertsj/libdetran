//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  WGSolver.hh
 *  @brief WGSolver class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_WGSOLVER_HH_
#define detran_WGSOLVER_HH_

#include "Solver.hh"
#include "angle/MomentToDiscrete.hh"
#include "transport/Sweeper.hh"
#include "transport/Sweeper2DMOC.hh"
#include "transport/SweepSource.hh"
#include "utilities/MathUtilities.hh"
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class WGSolver
 *  @brief Solve the within-group transport equation.
 *
 *  The within-group transport equation in operator form is
 *  @f[
 *      @mathbf{L}\psi =
 *        @mathbf{MS}\phi + Q
 *  @f]
 *  where \f$ \mathbf{L} \f$ is the streaming and collision operator,
 *  \mathbf{M} is the moment-to-discrete operator,
 *  \mathbf{S} is the scattering operator, and  Q represents any
 *  source considered fixed, which includes in-scatter, fission, and
 *  external sources.
 *
 *  What we are really after is the scalar flux and possibly its higher
 *  order moments.  Consequently, we are able to solve a somewhat different
 *  problem then the within group transport equation above.  Let us operate
 *  on both sides by \f$\mathbf{L}^{-1}\f$ followed
 *  by \f$ \mathbf{D}\f$ to get
 *  \f[
 *      (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi
 *        = \mathbf{D} \mathbf{L}^{-1} Q \, .
 *  \f]
 *  Here, \f$\mathbf{D}\f$ is the discrete-to-moment operator, defined
 *  such that \f$ \phi = \mathbf{D}\psi \f$.
 *
 *  Notice this is nothing but a linear system of the form
 *  \f$ \mathbf{A}x = b \f$ where
 *  \f[
 *      @mathbf{A} = (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})
 *  \f]
 *  and
 *  \f[
 *      b = \mathbf{D} \mathbf{L}^{-1} Q \, .
 *  \f]
 *  Moreover, \f$ b\f$ is just the uncollided flux.
 *
 *  A nice overview of approaches for this inner iteration is
 *  given by Larsen and Morel in <em> Nuclear Computational Science </em>.
 *
 *  Input parameters specific to WGSolver and derived classes:
 *    - inner_max_iters        [100]
 *    - inner_tolerance        [1e-5]
 *    - inner_print_out        [2], 0=never, 1=final, 2=every interval
 *    - inner_print_interval   [10]
 *
 *  @sa WGSolverSI, WGSolverGMRES
 */
//---------------------------------------------------------------------------//

template <class D>
class WGSolver: public Solver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<WGSolver>            SP_solver;
  typedef Solver<D>                                 Base;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_externalsource          SP_externalsource;
  typedef typename Base::vec_externalsource         vec_externalsource;
  typedef typename Base::SP_fissionsource 		    SP_fissionsource;
  typedef typename Sweeper<D>::SP_sweeper           SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource   SP_sweepsource;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef State::moments_type                       moments_type;
  typedef detran_angle::MomentToDiscrete::SP_MtoD   SP_MtoD;
  typedef typename Base::size_t                     size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  @param state             State vectors, etc.
   *  @param mat               Material definitions.
   *  @param quadrature        Angular mesh.
   *  @param boundary          Boundary fluxes.
   *  @param external_source   User-defined external source.
   *  @param fission_source    Fission source.
   */
  WGSolver(SP_state       	            state,
           SP_material                  material,
           SP_quadrature                quadrature,
           SP_boundary    	            boundary,
           const vec_externalsource    &q_e,
           SP_fissionsource             q_f,
           bool                         multiply);

  /// Virtual destructor
  virtual ~WGSolver(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL WITHIN-GROUP SOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// Solve the within group equation for group g
  virtual void solve(const size_t g) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Expose base members
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_boundary;
  using Base::d_externalsources;
  using Base::d_fissionsource;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;

  /// Angular mesh.
  SP_quadrature d_quadrature;
  /// Sweeper over the space-angle domain.
  SP_sweeper d_sweeper;
  /// Sweep source constructor
  SP_sweepsource d_sweepsource;
  /// Group we are solving.
  size_t d_g;
  /// Number of sweeps
  size_t d_number_sweeps;
  /// Flag for multiplying fixed source
  bool d_multiply;

private:

  /// Setup the templated sweeper
  bool set_sweep(std::string equation);

};

} // end namespace detran

#endif /* detran_WGSOLVER_HH_ */

//---------------------------------------------------------------------------//
//              end of WGSolver.hh
//---------------------------------------------------------------------------//
