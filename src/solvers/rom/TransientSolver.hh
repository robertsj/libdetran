/*
 * TransientSolver.hh
 *
 *  Created on: Sep 6, 2020
 *      Author: rabab
 */

#ifndef SOLVERS_ROM_TRANSIENTSOLVER_HH_
#define SOLVERS_ROM_TRANSIENTSOLVER_HH_

#include "FixedSourceManager.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/LinearSolver.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "kinetics/KineticsMaterial.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/vector/Vector.hh"
#include "callow/matrix/Matrix.hh"
#include "transport/State.hh"
#include "Kinetic_Mat.hh"
#include "kinetics/Precursors.hh"
#include "kinetics/SyntheticSource.hh"
#include "utilities/Definitions.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "OperatorProjection.hh"



using namespace detran;
using namespace callow;

class TransientSolver
{

public:
	typedef callow::MatrixDense::SP_matrix            SP_matrix;
	typedef callow::Vector::SP_vector                 SP_vector;
	typedef detran_utilities::InputDB::SP_input       SP_input;
	typedef detran_utilities::vec_int                 vec_int;
	typedef detran_utilities::vec_dbl                 vec_dbl;
    typedef detran_utilities::vec2_dbl                vec2_dbl;
	typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
    typedef TimeDependentMaterial::SP_material        SP_material;
	typedef Precursors::SP_precursors                 SP_precursors;
	typedef FissionSource::SP_fissionsource           SP_fissionsource;
	typedef State::SP_state                           SP_state;
    typedef SyntheticSource::vec_states               vec_states;
    typedef SyntheticSource::vec_precursors           vec_precursors;
    typedef DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
    typedef DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
    typedef LinearSolverCreator::SP_solver            SP_solver;

	TransientSolver(SP_input inp, SP_mesh mesh, SP_material material, SP_matrix flux_basis, SP_matrix precursors_basis);
	void initialize_precursors();
	void ProjectInitial(SP_state initial_state);
	void Solve(SP_state initial_state);
	void step(double t, double dt);
	void UpdatePrecursors(double dt);
	void UpdateOperator(double dt);
	void compute_b();

private:
  /// Vector of previous states
  vec_states d_states;
  vec_precursors d_precursors;

  /// Working precursor vector.
  SP_state d_state;
  SP_precursors d_precursor;

  /// Fission source
  SP_fissionsource d_fissionsource;

  SP_matrix d_F;
  SP_matrix d_P;
  SP_matrix d_G;
  SP_matrix d_D;
  SP_matrix d_L;
  SP_matrix s_A;
  SP_matrix d_B;

  /// Time step size
  double d_dt;

  // user settings
  SP_mesh d_mesh;
  SP_input d_inp;
  SP_material d_material;
  SP_matrix d_flux_basis;
  SP_matrix d_precursors_basis;
  SP_matrix d_A;

  // reduced vectors
  SP_vector d_P0;
  SP_vector d_phi0;
  SP_vector d_phi;

  SP_vector d_P_r;
  SP_vector d_P0_r;
  SP_vector d_phi0_r;
  SP_vector d_phi_r;
  SP_vector d_b;

  int d_num_cells;
  int d_number_groups;
  int d_precursors_group ;
  int d_rf;
  int d_rc;

  SP_solver d_solver;
};

#endif /* SOLVERS_ROM_TRANSIENTSOLVER_HH_ */
