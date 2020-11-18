//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TransientSolver.hh
 *  @brief TransientSolver class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef SOLVERS_ROM_TRANSIENTSOLVER_HH_
#define SOLVERS_ROM_TRANSIENTSOLVER_HH_

#include "callow/utils/Initialization.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/LinearSolver.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "kinetics/KineticsMaterial.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/vector/Vector.hh"
#include "transport/State.hh"
#include "kinetics/Precursors.hh"
#include "kinetics/SyntheticSource.hh"
#include "utilities/Definitions.hh"
#include "OperatorProjection.hh"
#include "KineticMatrices.hh"


using namespace detran;
using namespace callow;


class TransientSolver
{

public:
  typedef callow::MatrixDense::SP_matrix            SP_matrix;
  typedef callow::Vector::SP_vector                 SP_vector;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef TimeDependentMaterial::SP_material        SP_material;
  typedef Precursors::SP_precursors                 SP_precursors;
  typedef FissionSource::SP_fissionsource           SP_fissionsource;
  typedef State::SP_state                           SP_state;
  typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
  typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
  typedef LinearSolverCreator::SP_solver            SP_solver;
  typedef callow::Matrix::SP_matrix                 SP_matrix_Base;


  TransientSolver(SP_input inp, SP_mesh mesh, SP_material material, SP_matrix flux_basis, SP_matrix precursors_basis);

  void Solve(SP_state initial_state);

private:
  /// State vector
  SP_state d_state;
  /// Fission source
  SP_fissionsource d_fissionsource;
  /// Precursors production matrix
  SP_matrix d_precursors_production;
  /// Precursors decsy matrix
  SP_matrix d_precursors_decay;
  /// Delayed neutron production
  SP_matrix d_delayed_production;
  /// Gain matrix
  SP_matrix_Base d_G;
  /// Loss matrix
  SP_matrix_Base d_L;
  /// operator
  SP_matrix d_A;

  SP_matrix d_LF;
  SP_matrix d_Lr;
  SP_matrix d_Gr;
  /// Reduced solution
  SP_matrix d_sols;
  /// Flux solution matrix at all time steps
  SP_matrix d_flux;
  /// Precursors solution matrix at all time steps
  SP_matrix d_precursors;
  /// Reduced flux vector
  SP_matrix d_flux_r;
  /// Reduced precursors vector
  SP_matrix d_precursors_r;

  /// Time step size
  double d_dt;
  /// Final time
  double  d_final_time;
  /// Number of time steps
  double d_number_steps;

  /// Mesh
  SP_mesh d_mesh;
  /// Input
  SP_input d_inp;
  /// Solver setting
  SP_input db;
  /// Material
  SP_material d_material;
  /// Flux basis
  SP_matrix d_flux_basis;
  /// Precursors basis
  SP_matrix d_precursors_basis;
  /// Precursor vector in previous time step
  SP_vector d_P0;
  /// Flux vector in previous time step
  SP_vector d_phi0;
  /// Precursor vector in current time step
  SP_vector d_phi;
  /// Current vector of the reduced flux and precursors
  SP_vector d_sol_r;
  /// Vector of the reduced flux and precursors at the previous time step
  SP_vector d_sol0_r;
  /// The projected vector of the precursors
  SP_vector d_P_r;
  /// The projected vector of the initial precursors
  SP_vector d_P0_r;
  /// The projected vector of the initial flux
  SP_vector d_phi0_r;
  /// The projected vector of the flux
  SP_vector d_phi_r;
  ///
  SP_vector d_b;

  /// Number of cells
  int d_num_cells;
  /// number of energy groups
  int d_number_groups;
  /// Number of precursors group
  int d_precursors_group ;
  /// Flux rank
  int d_rf;
  /// Precursors rank
  int d_rc;
  /// Linear solver
  SP_solver d_solver;
  /// Compute the initial precursors concentration
  void initialize_precursors();
  /// Project the initial flux and precursors on space of the reduced basis
  void ProjectInitial();
  /// Construct matrix
  void Construct_Operator(double t, double dt);
  /// Update the operator 
  void Refersh_Operator();
  /// Reconstruct the full order solution
  void reconstruct();
};

#endif /* SOLVERS_ROM_TRANSIENTSOLVER_HH_ */
