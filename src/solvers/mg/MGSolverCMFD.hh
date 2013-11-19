//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverCMFD.hh
 *  @brief MGSolverCMFD class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGSolverCMFD_HH_
#define detran_MGSolverCMFD_HH_

#include "MGTransportSolver.hh"
#include "solvers/mg/CMFDLossOperator.hh"
#include "transport/CoarseMesh.hh"
#include "transport/CurrentTally.hh"
#include "transport/Homogenize.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class MGSolverCMFD
 *  @brief Solves multigroup problems via CMFD
 */
//----------------------------------------------------------------------------//
/**
 *  @example solvers/test/test_MGSolverCMFD.cc
 *
 *  Test of MGSolverCMFD class.
 */
template <class D>
class MGSolverCMFD: public MGTransportSolver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MGTransportSolver<D>                              Base;
  typedef typename Base::SP_solver                          SP_solver;
  typedef typename Base::SP_wg_solver                       SP_wg_solver;
  typedef typename Base::SP_input                           SP_input;
  typedef typename Base::SP_state                           SP_state;
  typedef typename Base::SP_mesh                            SP_mesh;
  typedef typename Base::SP_material                        SP_material;
  typedef typename Base::SP_quadrature                      SP_quadrature;
  typedef typename Base::SP_boundary                        SP_boundary;
  typedef typename Base::SP_externalsource                  SP_externalsource;
  typedef typename Base::vec_externalsource                 vec_externalsource;
  typedef typename Base::SP_fissionsource                   SP_fissionsource;
  typedef typename Base::size_t                             size_t;
  typedef detran_utilities::vec_dbl                         vec_dbl;
  typedef detran_utilities::vec2_dbl                        vec2_dbl;
  typedef detran_utilities::vec_size_t                      vec_size_t;
  typedef detran_utilities::vec_int                         vec_int;
  typedef CoarseMesh::SP_coarsemesh                         SP_coarsemesh;
  typedef typename CurrentTally<D>::SP_currenttally         SP_tally;
  typedef typename CMFDLossOperator<D>::SP_lossoperator     SP_operator;
  typedef callow::LinearSolver::SP_solver                   SP_linearsolver;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param q_e               Vector of user-defined external sources
   *  @param q_f               Fission source.
   *  @param multiply          Flag for a multiplying fixed source problem
   */
  MGSolverCMFD(SP_state                   state,
               SP_material                material,
               SP_boundary                boundary,
               const vec_externalsource  &q_e,
               SP_fissionsource           q_f,
               bool                       multiply = false);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the multigroup equations.
  void solve(const double keff = 1.0);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Return number of sweeps
  int number_sweeps() const
  {
    return d_wg_solver->get_sweeper()->number_sweeps();
  }

  /// Perform one sweep through all energies
  void energy_sweep();

  /// Getters
  SP_operator loss_operator() {return d_operator;}
  SP_mesh coarse_mesh() {return d_coarse_mesh;}
  SP_tally tally() {return d_tally;}
  SP_material coarse_material() {return d_coarse_material;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

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

  /// Lower and upper group bounds
  //@{
  int d_lower;
  int d_upper;
  //@}

  SP_coarsemesh d_coarsener;
  SP_mesh d_coarse_mesh;
  SP_tally d_tally;
  /// Linear solver for inverting the operator
  SP_linearsolver d_solver;
  /// Solver database
  SP_input d_solver_db;
  /// CMFD loss operator
  SP_operator d_operator;
  /// Diffusion coefficient homogenization option
  size_t d_diff_coef_weight;
  /// Relaxation factor
  double d_omega;
  /// Coarse material
  SP_material d_coarse_material;


  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  void update(const double keff);
  void compute_current();

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#endif /* detran_MGSolverCMFD_HH_ */

//----------------------------------------------------------------------------//
//              end of MGSolverCMFD.hh
//----------------------------------------------------------------------------//
