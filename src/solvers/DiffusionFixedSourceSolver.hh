//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DiffusionFixedSourceSolver.hh
 *  @brief  DiffusionFixedSourceSolver class definition
 *  @author Jeremy Roberts
 *  @date   Sep 11, 2012
 */
//---------------------------------------------------------------------------//

#ifndef DIFFUSIONFIXEDSOURCESOLVER_HH_
#define DIFFUSIONFIXEDSOURCESOLVER_HH_

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
 *  @class DiffusionFixedSourceSolver
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
class DiffusionFixedSourceSolver
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum fixed_types
  {
    FIXED, MULTIPLY, ITERATE, END_FIXED_TYPES
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  //
  typedef BoundaryDiffusion<D>                        Boundary_T;
  typedef typename Boundary_T::SP_boundary            SP_boundary;
  //
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_source;
  //
  typedef State::moments_type                         moments_type;
  //
  typedef DiffusionLossOperator::SP_lossoperator      SP_lossoperator;
  typedef DiffusionGainOperator::SP_gainoperator      SP_gainoperator;
  typedef callow::Vector                              Vector_T;
  typedef Vector_T::SP_vector                         SP_vector;
  typedef callow::LinearSolverCreator                 Creator_T;
  typedef callow::LinearSolver::SP_solver             SP_solver;
  typedef callow::LinearSolver::SP_preconditioner     SP_preconditioner;
  typedef callow::MatrixBase::SP_matrix               SP_matrix;
  //
  typedef detran_utilities::size_t                    size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   *  @param state      state vector
   */
  DiffusionFixedSourceSolver(SP_input input,
                             SP_material material,
                             SP_mesh mesh,
                             SP_state state);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Build the right hand side
   *
   *  The source consists of external volume and boundary sources.
   *  Currently, the client must set the boundary current directly for
   *  a boundary source.
   *
   *  Note, the volume source is optional.  If it is not present, only
   *  the boundary current contributes.
   *
   *  @param q  Pointer to volume source
   */
  void build_source(SP_source q = SP_source(0));

  /**
   *  @brief Set the right hand side directly
   *
   *  This allows the client to set an explicit right hand side
   *  for the solve, useful when the diffusion solver is
   *  nested within another iteration.  Note that no boundary
   *  source is constructed from the current boundary information.
   *
   *  @param b  Pointer to right hand side vector
   */
  void build_source(SP_vector b);

  /// Solve the fixed source diffusion problem
  void solve();

  /// @name Getters
  /// @{
  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }
  SP_matrix lossoperator() const { return d_M; }
  SP_matrix gainoperator() const { return d_F; }
  SP_vector phi() { return d_phi; }
  SP_vector Q() { return d_Q; }
  SP_vector Q_total() { return d_Q_total; }
  /// @}

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// User input
  SP_input d_input;
  /// Material database
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// State vector
  SP_state d_state;
  /// Boundary container
  SP_boundary d_boundary;
  /// Loss operator
  SP_lossoperator d_M;
  /// Fission operator
  SP_gainoperator d_F;
  /// Preconditioner
  SP_preconditioner d_P;
  /// Unknown vector
  SP_vector d_phi;
  /// Saved unknown vector
  SP_vector d_phi_old;
  /// Fixed source
  SP_vector d_Q;
  /// Fixed source plus fission
  SP_vector d_Q_total;
  /// Number of unknowns
  size_t d_problem_size;
  /// Solver type (e.g. gmres)
  std::string d_solver_type;
  /// Type indicating how to treat fission
  size_t d_fixed_type;
  /// Maximum linear solver iterations
  size_t d_maximum_iterations;
  /// Maximum fission iterations
  size_t d_maximum_fission_iterations;
  /// Tolerance for linear solve
  double d_tolerance;
  /// Tolerance for the fission iteration solve
  double d_fission_tolerance;
  /// Scaling factor for the gains operator (i.e. 1/keff)
  double d_fission_scaling;
  /// Pointer to linear solver
  SP_solver d_solver;
  int d_print_out;
  int d_print_interval;
  bool d_adjoint;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Solve directly via M*phi = Q
  void solve_fixed();
  /// Solve with iteration on the fission source, M*phi(n+1) = (F/k)*phi(n) + Q
  void solve_iterate();

  /// Build the volume source
  void build_volume_source(SP_source q);
  /// Build the boundary source
  void build_boundary_source();

  /// Fill the boundary with outgoing current
  void fill_boundary();

};

} // end namespace detran

#endif // DIFFUSIONFIXEDSOURCESOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file DiffusionFixedSourceSolver.hh
//---------------------------------------------------------------------------//
