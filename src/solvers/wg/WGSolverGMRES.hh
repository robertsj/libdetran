//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  WGSolverGMRES.hh
 *  @brief WGSolverGMRES class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_WGSOLVERGMRES_HH_
#define detran_WGSOLVERGMRES_HH_

#include "WGSolver.hh"
#include "WGTransportOperator.hh"
#include "WGPreconditioner.hh"
#include "callow/solver/LinearSolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class WGSolverGMRES
 *  @brief Solve the within-group problem with GMRES.
 *
 *  From @ref WGSolver, we know the within-group problem can be
 *  written in operator notation as
 *  @f[
 *      (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi
 *      = \mathbf{D} \mathbf{L}^{-1} Q \, ,
 *  @f]
 *  or
 *  @f[
 *      \mathbf{A}x = b \, .
 *  @f]
 *
 *  This class couples with callow to make available its set of applicable
 *  solvers, the default being GMRES.  Other solvers are selected
 *  by the parameter database for callow solvers.  PETSc solvers are available
 *  through the callow interface as well.
 *
 *  For Krylov iterations to perform successfully, preconditioning
 *  is often required.  A good preconditioner @f$ \mathbf{M} @f$
 *  is in some way "similar" to the operator @f$ \mathbf{A} @f$, and
 *  applying its inverse @f$ \mathbf{M}^{-1} @f$ can be done cheaply.
 */
//---------------------------------------------------------------------------//

template <class D>
class WGSolverGMRES: public WGSolver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

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
  typedef detran_utilities::vec_dbl             vec_dbl;
  //
  typedef WGTransportOperator<D>                Operator_T;
  typedef typename Operator_T::SP_operator      SP_operator;
  typedef WGPreconditioner::SP_preconditioner   SP_preconditioner;
  typedef callow::LinearSolver::SP_solver       SP_linearsolver;
  typedef callow::Vector::SP_vector             SP_vector;

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
  WGSolverGMRES(SP_state                   state,
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
  using Base::d_g;

  /// Main linear solver
  SP_linearsolver d_solver;
  /// Preconditioner
  SP_preconditioner d_pc;
  /// Operator "A" in "Ax = b"
  SP_operator d_operator;
  /// Solution vector
  SP_vector d_x;
  /// Right hand side
  SP_vector d_b;
  /// Number of reflective iterations
  int d_reflective_solve_iterations;
  /// Flag to update the angular fluxes, including cell and boundary
  bool d_update_angular_flux;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Build the right hand side.
  void build_rhs(State::moments_type &B);

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "WGSolverGMRES.i.hh"

#endif /* detran_WGSOLVERGMRES_HH_ */

//----------------------------------------------------------------------------//
//              end of WGSolverGMRES.hh
//----------------------------------------------------------------------------//
