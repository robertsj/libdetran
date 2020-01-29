//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenDiffusion.hh
 *  @brief EigenDiffusion class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENDIFFUSION_HH_
#define detran_EIGENDIFFUSION_HH_

#include "Eigensolver.hh"
#include "DiffusionGainOperator.hh"
#include "DiffusionLossOperator.hh"
#include "callow/solver/EigenSolverCreator.hh"

namespace detran
{


/**
 *  @class EigenDiffusion
 *  @brief Solves the diffusion eigenvalue directly
 *
 *  While the other eigensolvers can be used in conjuction with the
 *  multigroup diffusion solver, in this eigensolver, both the loss
 *  and gains operator are explicitly constructed for solution via
 *  a generalized eigensolver.  Including this capability is mostly
 *  for benchmarking the performance of the other approach.
 *
 */
template <class D>
class EigenDiffusion: public Eigensolver<D>
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
  typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  typedef callow::EigenSolverCreator                Creator_T;
  typedef callow::EigenSolver::SP_solver            SP_eigensolver;
  typedef callow::Vector::SP_vector                 SP_vector;
  typedef detran_utilities::size_t                  size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver         Multigroup solver
   */
  EigenDiffusion(SP_mg_solver mg_solver);

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

  /// Gain operator
  SP_gainoperator d_F;
  /// Loss operator
  SP_lossoperator d_M;
  /// Flux vector
  SP_vector d_phi;
  /// Working vector
  SP_vector d_work;
  /// Eigenvalue solver
  SP_eigensolver d_eigensolver;

};

} // namespace detran

#endif // detran_EIGENDIFFUSION_HH_

//----------------------------------------------------------------------------//
//              end of file EigenDiffusion.hh
//----------------------------------------------------------------------------//
