//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   WGSolverSI.hh
 *  @brief  WGSolverSI class definition.
 *  @note   Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_WGSOLVERSI_HH_
#define detran_WGSOLVERSI_HH_

// Detran
#include "WGSolver.hh"

#include <iostream>

namespace detran
{

/**
 *  @class WGSolverSI
 *  @brief Solve the within-group transport equation via source iteration.
 *
 *  Source iteration (or Richardson iteration) solves the transport
 *  equation (see \ref WGSolver)
 *
 *  \f[
 *      (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi
 *        = \mathbf{D} \mathbf{L}^{-1} Q \, .
 *  \f]
 *
 *  via the process
 *
 *  \f[
 *      \phi^{(n+1)} = \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi^{(n)}
 *        + \mathbf{D} \mathbf{L}^{-1} Q \, .
 *  \f]
 *
 *  This is a "hand-coded" implementation.  Essentially the same solver
 *  can be had via callow's Richardson iteration (and via callow's
 *  interface to PETSc's Richardson).
 *
 *  This implementation is useful because it provides access to various
 *  nonlinear acceleration schemes not applicable to nonstationary solvers.
 */

template <class D>
class WGSolverSI: public WGSolver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef WGSolver<D>                           Base;
  typedef typename Base::SP_solver              SP_solver;
  typedef typename Base::SP_input               SP_input;
  typedef typename Base::SP_state               SP_state;
  typedef typename Base::SP_mesh                SP_mesh;
  typedef typename Base::SP_material            SP_material;
  typedef typename Base::SP_quadrature          SP_quadrature;
  typedef typename Base::SP_boundary            SP_boundary;
  typedef typename Base::SP_MtoD                SP_MtoD;
  typedef typename Base::SP_externalsource      SP_externalsource;
  typedef typename Base::vec_externalsource     vec_externalsource;
  typedef typename Base::SP_fissionsource       SP_fissionsource;
  typedef typename Base::SP_sweeper             SP_sweeper;
  typedef typename Base::SP_sweepsource         SP_sweepsource;
  typedef typename Base::moments_type           moments_type;
  typedef typename Base::size_t                 size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param mat               Material definitions.
   *  @param quadrature        Angular mesh.
   *  @param boundary          Boundary fluxes.
   *  @param external_source   User-defined external source.
   *  @param fission_source    Fission source.
   *  @param multiply          Flag for fixed source multiplying problem
   */
  WGSolverSI(SP_state                   state,
             SP_material                material,
             SP_quadrature              quadrature,
             SP_boundary                boundary,
             const vec_externalsource  &q_e,
             SP_fissionsource           q_f,
             bool                       multiply);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL WITHIN-GROUP SOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the within group equation.
  void solve(const size_t g);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  // Make inherited data visible
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_quadrature;
  using Base::d_boundary;
  using Base::d_sweeper;
  using Base::d_sweepsource;
  using Base::d_tolerance;
  using Base::d_maximum_iterations;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;
  using Base::d_g;

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "WGSolverSI.i.hh"

#endif /* detran_WGSOLVERSI_HH_ */

//----------------------------------------------------------------------------//
//              end of WGSolverSI.hh
//----------------------------------------------------------------------------//
