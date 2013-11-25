//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverGMRES.hh
 *  @brief MGSolverGMRES class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGSOLVERGMRES_HH_
#define detran_MGSOLVERGMRES_HH_

#include "MGTransportSolver.hh"
#include "MGTransportOperator.hh"
#include "MGPreconditioner.hh"
#include "callow/solver/LinearSolver.hh"
#include "transport/Sweeper.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class MGSolverGMRES
 *  @brief Solves the multigroup transport equation via GMRES.
 *
 *  Traditionally, the Gauss-Seidel method has been used for multigroup
 *  problems. For each group, the within-group equation is solved, and
 *  the the fluxes are updated for use in the next group.  However, for
 *  problems with significant upscatter, Gauss-Seidel can be quite
 *  expensive, even when GMRES (or some better-than-source-iteration
 *  scheme) is used for the within group solve.  As an
 *  alternative, we can apply GMRES (or other Krylov solvers) to
 *  the multigroup
 *  problem directly.  The linear system is then
 *  \f[
 *     \left( ( \mathbf{I} -
 *         \left(\begin{array}{ccc}
 *             T_1  & \cdots & 0     \\
 *             0    & \ddots & 0     \\
 *             0    & 0      & T_G
 *         \end{array}\right) \cdot
 *         \left(\begin{array}{ccc}
 *             M    & \cdots & 0     \\
 *             0    & \ddots & 0     \\
 *             0    & 0      & M
 *         \end{array}\right) \cdot
 *         \left(\begin{array}{ccc}
 *            \mathbf{S}_{11} & \cdots & \mathbf{S}_{1G} \\
 *            \vdots          & \ddots & \vdots          \\
 *            \mathbf{S}_{G1} & 0      & \mathbf{S}_{GG}
 *         \end{array}\right)
 *     \right ) \cdot
 *     \left[ \begin{array}{c}
 *         \phi_1 \\
 *         \vdots \\
 *         \phi_G
 *     \end{array} \right] =
 *     \left[ \begin{array}{c}
 *          \mathbf{T}_1 q_1 \\
 *          \vdots \\
 *          \mathbf{T}_G q_G
 *     \end{array} \right] \, .
 *  \f]
 *  Of course, this can be written succinctly in the same way as the within-
 *  group equation:
 *  \f[
 *     (\mathbf{I}-\mathbf{TMS})\phi = \mathbf{T}q \, ,
 *  \f]
 *  where \f$ \mathbf{T} = D\mathbf{L}^{-1} \f$ is the sweeping operator with
 *  moment contributions added implicitly, and where the Krylov vectors are
 *  energy-dependent.
 *
 *  By default, only the energy block in which upscatter occurs is
 *  solved via Krylov methods.  Because Gauss-Seidel is exact for
 *  downscatter, it is used for the downscatter-only block.  The user can
 *  switch this using "outer_upscatter_cutoff".
 *
 *  Reference:
 *    Evans, T., Davidson, G. and Mosher, S. "Parallel Algorithms for
 *    Fixed-Source and Eigenvalue Problems", NSTD Seminar (ORNL), May 27, 2010.
 *
 *  @note While the intent of this class is to provide GMRES as a multigroup
 *        solver, it makes available any of the applicable \ref callow
 *        solvers, included PETSc solvers when enabled.  The default solver
 *        is the built-in GMRES, which is not a highly-tuned implementation!
 *
 *  \todo Consider better ways to handle memory between Vec and vector.  May
 *        want to devise a moment container based on pointer that allows
 *        one to swap memory temporarily (as done with Vec)
 *
 */
//----------------------------------------------------------------------------//
/**
 *  @example solvers/test/test_MGSolverGMRES.cc
 *
 *  Test of MGSolverGMRES class.
 */
template <class D>
class MGSolverGMRES: public MGTransportSolver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGTransportSolver<D>                      Base;
  typedef typename Base::SP_solver                  SP_solver;
  typedef typename Base::SP_wg_solver               SP_wg_solver;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_quadrature              SP_quadrature;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_externalsource          SP_externalsource;
  typedef typename Base::vec_externalsource         vec_externalsource;
  typedef typename Base::SP_fissionsource           SP_fissionsource;
  typedef typename Sweeper<D>::SP_sweeper           SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource   SP_sweepsource;
  typedef callow::LinearSolver::SP_solver           SP_linearsolver;
  typedef MGTransportOperator<D>                    Operator_T;
  typedef typename Operator_T::SP_operator          SP_operator;
  typedef callow::Vector::SP_vector                 SP_vector;
  typedef MGPreconditioner::SP_preconditioner       SP_pc;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::size_t                  size_t;
  typedef detran_utilities::vec_size_t              groups_t;
  typedef groups_t::iterator                        groups_iter;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param q_e               Vector of user-defined external sources
   *  @param q_f               Fission source.
   *  @param multiply          Flag for a multiplying fixed source problem
   */
  MGSolverGMRES(SP_state                   state,
                SP_material                material,
                SP_boundary                boundary,
                const vec_externalsource  &q_e,
                SP_fissionsource           q_f,
                bool                       multiply = false);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// Solve the multigroup equations.
  void solve(const double keff = 1.0);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Return number of sweeps
  int number_sweeps() const;

  /// Return the transport operator
  SP_operator get_operator();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Expose base members.
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_quadrature;
  using Base::d_boundary;
  using Base::d_externalsources;
  using Base::d_fissionsource;
  using Base::d_downscatter;
  using Base::d_number_groups;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;
  using Base::d_wg_solver;
  using Base::d_multiply;

  /// Main linear solver
  SP_linearsolver d_solver;
  /// Preconditioner
  SP_pc d_pc;
  /// Operator "A" in "Ax = b"
  SP_operator d_operator;
  /// Solution vector
  SP_vector d_x;
  /// Right hand side
  SP_vector d_b;
  /// Size of the moments portion of d_X
  int d_moments_size;
  /// Size of the moments portion of d_X in a group
  int d_moments_size_group;
  /// Size of the boundary portion of d_X
  int d_boundary_size;
  /// Size of the boundary portion of d_X in a group
  int d_boundary_size_group;

  /**
   *  @brief Only groups equal to or above this cutoff are
   *         subject to Krylov iterations.
   *
   *  While \ref Material has an upscatter cutoff that it computes
   *  internally based on the data, the user can set this cutoff
   *  to a different value for solving.  By default, the solver
   *  cutoff is equal to the \ref Material cutoff, and so those
   *  groups into which no upscatter occurs are solved by
   *  Gauss-Seidel, and the remaining groups are solved by the
   *  PETSc solver selected.  The user can set the cutoff to
   *  zero to use the PETSc solver for the entire energy spectrum,
   *  or any value in between zero through the actual cutoff.
   *  The user \e cannot set the cutoff any higher than the \ref
   *  Material cutoff, since that would be a different problem.
   *
   */
  int d_krylov_group_cutoff;
  int d_lower;
  int d_upper;
  /// Krylov block size (number of groups in Krylov portion of solve)
  int d_number_active_groups;
  /// Count of reflective solve iterations
  int d_reflective_solve_iterations;
  /// Sweeper
  SP_sweeper d_sweeper;
  /// Sweep source
  SP_sweepsource d_sweepsource;
  /// Flag to update the angular flux, including boundary and cell values
  bool d_update_angular_flux;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Build the right hand side.
  void build_rhs(State::vec_moments_type &B);

  /**
   *  @brief Sweep through the upscatter block once
   *
   *  For computing the right hand side (i.e. the uncollided flux), we
   *  need to sweep through all the groups once.  Within a group, we need
   *  to iterate on reflecting conditions.  This process is also need at
   *  the end of a solve to pick up the boundary fluxes.
   *
   *  @param phi    multigroup fluxes to update
   */
  void group_sweep(State::vec_moments_type &phi);

  /// Solve the within-group reflection problem
  void solve_wg_reflection(size_t g, State::moments_type &phi_g);

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "MGSolverGMRES.i.hh"

#endif /* detran_MGSOLVERGMRES_HH_ */

//----------------------------------------------------------------------------//
//              end of MGSolverGMRES.hh
//----------------------------------------------------------------------------//
