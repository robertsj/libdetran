//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenGD.hh
 *  @brief EigenGD class definition.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENGD_HH_
#define detran_EIGENGD_HH_

#include "Eigensolver.hh"
#include "EnergyDependentEigenLHS.hh"
#include "solvers/mg/MGTransportOperator.hh"
#include "callow/solver/EigenSolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class EigenGD
 *  @brief Solves the eigenvalue problem using a generalized Davidson method
 *
 *  The multigroup eigenvalue problem can be written compactly as
 *  @f[
 *     \mathbf{A}\phi = \frac{1}{k} \mathbf{F} \phi \,
 *  @f]
 *  following the discussion in @ref Eigensolver.
 *
 *  By solving the generalized eigenvalue problem, we can avoid explicitly
 *  solving the multigroup transport problem and the sweeps it entails.
 */
//----------------------------------------------------------------------------//

template <class D>
class EigenGD: public Eigensolver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Eigensolver<D>                            Base;
  typedef typename Base::SP_solver                  SP_solver;
  typedef typename Base::SP_mg_solver               SP_mg_solver;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_fissionsource           SP_fissionsource;
  typedef callow::EigenSolver::SP_solver            SP_eigensolver;
  typedef MGTransportOperator<D>                    RHS_Operator_T;
  typedef EnergyDependentEigenLHS<D>                LHS_Operator_T;
  typedef callow::Vector::SP_vector                 SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver         Multigroup solver
   */
  EigenGD(SP_mg_solver mg_solver);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the eigenvalue problem.
  void solve();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  void update_preconditioner();

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Expose base members
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_fissionsource;
  using Base::d_number_groups;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;
  using Base::d_mg_solver;

  /// Eigensolver here callow's GD
  SP_eigensolver d_eigensolver;
  /// Operators
  //@{
  typename RHS_Operator_T::SP_operator d_A;
  typename LHS_Operator_T::SP_operator d_F;
  //@}
  /// Solution vector
  SP_vector d_x;
  /// Initial guess
  SP_vector d_x0;

};

} // namespace detran

#endif /* detran_EIGENGD_HH_ */

//----------------------------------------------------------------------------//
//              end of EigenGD.cc
//----------------------------------------------------------------------------//
