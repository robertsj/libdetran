//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenArnoldi.hh
 *  @brief EigenArnoldi class definition.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENSLEPC_HH_
#define detran_EIGENSLEPC_HH_

#include "Eigensolver.hh"
#include "EnergyIndependentEigenOperator.hh"
#include "callow/solver/EigenSolver.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class EigenArnoldi
 *  @brief Solves the eigenvalue problem using the Arnoldi method
 *
 *  The eigenvalue problem can be cast in the form
 *  @f[
 *     \mathbf{A}d = kd \,
 *  @f]
 *  where \f$ d \f$ is the fission density and \f$ k \f$ is the
 *  eigenvalue.  See Eigensolver for more details on this formulation.
 *
 *  As mentioned in @ref Eigensolver, Krylov methods require an
 *  adequately-converged multigroup solve.  The paper by Warsa et al.
 *  talks more about this, and it seems the multigroup
 *  convergence is most
 *  important in the first several iterations.  Thus, it might
 *  be worth implementing a dynamic tolerance at some point.
 *
 *  @note Detran must be configured with SLEPc for Arnoldi to be
 *        available via the callow interface.  Otherwise, the
 *        callow implementation of power iteration is used.
 *
 */
//----------------------------------------------------------------------------//

template <class D>
class EigenArnoldi: public Eigensolver<D>
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
  typedef EnergyIndependentEigenOperator<D>         Operator_T;
  typedef typename Operator_T::SP_operator          SP_operator;
  typedef callow::Vector::SP_vector                 SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver         Multigroup solver
   */
  EigenArnoldi(SP_mg_solver mg_solver);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the eigenvalue problem.
  void solve();

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

  /// Main linear solver
  SP_eigensolver d_eigensolver;
  /// Operator "A" in "Ax = b"
  SP_operator d_operator;
  /// Solution vector
  SP_vector d_x;
  /// Initial guess
  SP_vector d_x0;

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "EigenArnoldi.i.hh"

#endif /* detran_EIGENARNOLDI_HH_ */

//----------------------------------------------------------------------------//
//              end of EigenArnoldi.cc
//----------------------------------------------------------------------------//
