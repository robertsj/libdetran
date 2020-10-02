//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ROMSolver.hh
 *  @brief RoMSolver class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef SOLVERS_ROM_ROMSOLVER_HH_
#define SOLVERS_ROM_ROMSOLVER_HH_


#include "callow/vector/Vector.hh"
#include "callow/matrix/MatrixBase.hh"
#include "utilities/InputDB.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "FixedSourceManager.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/mg/DiffusionLossOperator.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "solvers/eigen/EnergyDependentEigenLHS.hh"
#include "solvers/eigen/EnergyIndependentEigenOperator.hh"
#include "solvers/mg/MGTransportOperator.hh"
#include "solvers/mg/MGSolverGMRES.hh"
#include "EigenGD.hh"
#include "solvers/mg/MGSolverGMRES.hh"
#include "OperatorProjection.hh"
#include <string>
#include <math.h>

using namespace detran;


//----------------------------------------------------------------------------//
/**
 *  @class ROMSolver
 *  @brief Solves the reduced eigenvalue problem
 */
//----------------------------------------------------------------------------//


template <class D>
class ROMSolver
{
public:
   typedef callow::EigenSolver::SP_solver            SP_eigensolver;
   typedef callow::EigenSolverCreator                Creator_T;
   typedef callow::MatrixBase::SP_matrix             SP_matrix;
   typedef callow::Vector::SP_vector                 SP_vector;
   typedef detran_utilities::InputDB::SP_input       SP_input;
   typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
   typedef detran_material::Material::SP_material        SP_material;
   typedef detran::DiffusionLossOperator::SP_lossoperator    SP_lossoperator;
   typedef detran::DiffusionGainOperator::SP_gainoperator    SP_gainoperator;
   typedef typename Eigensolver<D>::Fixed_T              Fixed_T;
   typedef typename Eigensolver<D>::SP_solver            SP_solver;
   typedef typename Fixed_T::SP_manager                  SP_mg_solver;
   typedef EnergyDependentEigenLHS<D>                LHS_Operator_T;
   typedef MGTransportOperator<D>                    RHS_Operator_T;
   typedef EnergyIndependentEigenOperator<D>             Operator_T;

   /// Constructor
   ROMSolver(SP_input inp, SP_mesh Mesh, SP_material mat);

   /// Solve the eigenvalue problem
   void Solve(SP_matrix d_U, SP_vector sol);

   /// Give the fundamental eigenvalue
   double keff()
   {
     return d_keff;
   }


private:
  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Input
  SP_input d_input;
  /// Mesh
  SP_mesh d_mesh;
  /// Material
  SP_material d_mat;
  /// Operator "A" in "Ax = \lambda B"
  SP_matrix  d_A;
  /// Operator "B" in "Ax = \lambda B"
  SP_matrix d_B;
  /// Basis rank
  int d_r;
  /// Operator size
  int d_n;
  /// Eigenvalue
  double d_keff;
  /// Eigen-vector of the reduced eigenvalue problem
  SP_vector x_rom;
  /// Reconstructed eigen-vector
  SP_vector x_fom;
  /// name of the problem operator
  std::string d_operator;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Sets the basis
  void SetBasis();
  /// Sets the operators of the problem
  void Set_FullOperators();
};

#endif /* SOLVERS_ROM_ROMSOLVER_HH_ */
