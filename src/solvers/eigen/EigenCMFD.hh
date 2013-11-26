//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenCMFD.hh
 *  @brief EigenCMFD class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENCMFD_HH_
#define detran_EIGENCMFD_HH_

#include "Eigensolver.hh"
#include "solvers/mg/MGSolverCMFD.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "callow/solver/EigenSolver.hh"
#include "callow/matrix/MatrixBase.hh"
#include "utilities/Definitions.hh"

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class EigenCMFD
 *  @brief Solves the eigenvalue problem via CMFD-accelerated power iteration
 */
//----------------------------------------------------------------------------//

template <class D>
class EigenCMFD: public Eigensolver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Eigensolver<D>                                Base;
  typedef typename Base::SP_solver                      SP_solver;
  typedef typename Base::SP_mg_solver                   SP_mg_solver;
  typedef typename Base::SP_input                       SP_input;
  typedef typename Base::SP_state                       SP_state;
  typedef typename Base::SP_mesh                        SP_mesh;
  typedef typename Base::SP_material                    SP_material;
  typedef typename Base::SP_boundary                    SP_boundary;
  typedef typename Base::SP_fissionsource               SP_fissionsource;
  typedef callow::EigenSolver::SP_solver                SP_eigensolver;
  typedef DiffusionGainOperator::SP_gainoperator        SP_gainoperator;
  typedef typename CMFDLossOperator<D>::SP_lossoperator SP_lossoperator;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver         Multigroup solver
   */
  EigenCMFD(SP_mg_solver mg_solver);

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

  /// Diffusion eigenvalue solver
  SP_eigensolver d_eigensolver;
  /// Eigensolver parameter database
  SP_input d_eigen_db;
  /// Over-relaxation parameter
  double d_omega;
  /// Gain operator
  SP_gainoperator d_gain;
  /// Loss operator
  SP_lossoperator d_loss;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  double cmfd_update();
  void compute_current();

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

//#include "EigenCMFD.i.hh"

#endif /* detran_EIGENCMFD_HH_ */

//----------------------------------------------------------------------------//
//              end of file EigenCMFD.hh
//----------------------------------------------------------------------------//
