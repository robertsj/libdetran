//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ReflectiveSolver.hh
 * \brief  ReflectiveSolver 
 * \author Jeremy Roberts
 * \date   Nov 2, 2012
 */
//---------------------------------------------------------------------------//

#ifndef REFLECTIVESOLVER_HH_
#define REFLECTIVESOLVER_HH_

#include "Solver.hh"
#include "SweepOperator.hh"
#include "callow/solver/LinearSolver.hh"

namespace detran
{


//---------------------------------------------------------------------------//
/**
 *  @class ReflectiveSolver
 *  @brief Solve the reflective boundary condition problem
 *
 *  This solves
 *  @f[
 *     \mathbf{L}\psi = q \, ,
 *  @f]
 *  with @f$ \psi @f$ subject to reflective conditions.  In all but
 *  pure reflection, this should take very few sweeps.  For pure
 *  vacuum conditions, just one sweep is required.
 *
 */
//---------------------------------------------------------------------------//

template <class D>
class ReflectiveSolver: public Solver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Solver<D>                             Base;
  typedef typename Base::SP_solver              SP_solver;
  typedef typename Base::SP_input               SP_input;
  typedef typename Base::SP_state               SP_state;
  typedef typename Base::SP_mesh                SP_mesh;
  typedef typename Base::SP_boundary            SP_boundary;
  typedef typename Base::SP_sweeper             SP_sweeper;
  typedef typename Base::SP_sweepsource         SP_sweepsource;
  typedef typename Base::moments_type           moments_type;
  typedef typename Base::size_t                 size_t;
  typedef detran_utilities::vec_dbl             vec_dbl;
  typedef callow::LinearSolver::SP_solver       SP_linearsolver;
  typedef SweepOperator<D>                      Operator_T;
  typedef typename Operator_T::SP_operator      SP_operator;
  typedef callow::Vector::SP_vector             SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param boundary          Boundary fluxes.
   *  @param sweeper
   *  @param source
   */
  ReflectiveSolver(SP_state        state,
                   SP_boundary     boundary,
                   SP_sweeper      sweeper,
                   SP_sweepsource  source);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Solve the reflection equation given flux moments
  void solve(SP_vector phi);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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
  /// Operator "A" in "Ax = b"
  SP_operator d_operator;
  /// Solution vector
  SP_vector d_x;
  /// Right hand side
  SP_vector d_b;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Set the templated operator function.
  PetscErrorCode set_operation();

  /// Build the right hand side.
  void build_rhs(State::moments_type &B);

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//#include "ReflectiveSolver.i.hh"

#endif // REFLECTIVESOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file ReflectiveSolver.hh
//---------------------------------------------------------------------------//
