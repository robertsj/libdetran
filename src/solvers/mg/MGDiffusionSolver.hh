//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGDiffusionSolver.hh
 *  @brief MGDiffusionSolver class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGDIFFUSIONSOLVER_HH_
#define detran_MGDIFFUSIONSOLVER_HH_

#include "MGSolver.hh"
#include "DiffusionLossOperator.hh"
#include "DiffusionGainOperator.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "external_source/ExternalSource.hh"
#include "transport/State.hh"
#include "utilities/Definitions.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/solver/LinearSolverCreator.hh"

namespace detran
{

/**
 *  @class MGDiffusionSolver
 *  @brief Solve a fixed source problem using diffusion
 *
 *  Fixed source problems can be solved in the following modalities:
 *    - fixed, no multiplication
 *    - fixed, multiplying via implicit fission (i.e. in the loss matrix)
 *    - fixed, multiplying via fission iteration
 *
 *  Mathematically, consider the diffusion equation in
 *  operator form:
 *  @f[
 *      T \phi = \frac{1}{k} F \phi + Q \, .
 *  @f]
 *  The first case assumes \f$ F = 0 \f$.  The second case
 *  uses the modified operator \f$ T' = T -  \frac{1}{k} F \phi \f$.
 *  The third case performs the iteration
 *  @f[
 *      T \phi^{n+1} = \frac{1}{k}F\phi^{n} + Q \, .
 *  @f]
 *  This third case is less efficient than the second case, but it
 *  allows one to pull out the solution following each fission
 *  iteration.  The number of such fission iterations can be limited
 *  by the user.
 *
 *  These three cases are selected via diffusion_fixed_type 0,1,2
 */
template <class D>
class MGDiffusionSolver: public MGSolver<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  // From base
  typedef MGSolver<D>                               Base;
  typedef typename Base::SP_solver                  SP_solver;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_externalsource          SP_externalsource;
  typedef typename Base::vec_externalsource         vec_externalsource;
  typedef typename Base::SP_fissionsource           SP_fissionsource;
  typedef State::moments_type                       moments_type;
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef callow::Vector                            Vector_T;
  typedef Vector_T::SP_vector                       SP_vector;
  typedef callow::LinearSolverCreator               Creator_T;
  typedef callow::LinearSolver::SP_solver           SP_linearsolver;
  typedef callow::LinearSolver::SP_preconditioner   SP_preconditioner;
  typedef callow::MatrixBase::SP_matrix             SP_matrix;
  typedef detran_utilities::size_t                  size_t;
  typedef detran_utilities::vec_int                 vec_int;
  typedef BoundaryDiffusion<D>                      Boundary_T;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param external_source   User-defined external source.
   *  @param fission_source    Fission source.
   *  @param multiply          Flag for multiplying fixed source problem
   */
  MGDiffusionSolver(SP_state                  state,
                    SP_material               material,
                    SP_boundary               boundary,
                    const vec_externalsource &q_e,
                    SP_fissionsource          q_f,
                    bool                      multiply);

  /// Refresh the solver.
  void refresh();

  /// Return the lossoperator
  SP_lossoperator lossoperator()
  {
    return d_M;
  }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve the fixed source diffusion problem
  void solve(const double keff = 1.0);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Expose base members
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
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
  using Base::d_multiply;

  /// Loss operator
  SP_lossoperator d_M;
  /// Fission operator
  SP_gainoperator d_F;
  /// Preconditioner
  SP_preconditioner d_P;
  /// Unknown vector
  SP_vector d_phi;
  /// Fixed source
  SP_vector d_Q;
  /// Number of unknowns
  size_t d_problem_size;
  /// Solver type (e.g. gmres)
  std::string d_solver_type;
  /// Pointer to linear solver
  SP_linearsolver d_solver;
  /// Scaling factor for fission source
  double d_keff;
  /// Boundary fill flag
  bool d_fill_boundary;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Build the volume source from external sources and fission
  void build_volume_source();
  /// Build the boundary source from the latest boundary values
  void build_boundary_source();
  /// Fill the state with the flux
  void fill_state();
  /// Fill the state with the current
  void fill_current();
  /// Fill the boundary with outgoing current
  void fill_boundary();

};

} // end namespace detran

#endif // detran_MGDIFFUSIONSOLVER_HH_

//----------------------------------------------------------------------------//
//              end of file MGDiffusionSolver.hh
//----------------------------------------------------------------------------//
