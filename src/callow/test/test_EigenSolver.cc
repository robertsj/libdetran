//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_EigenSolver.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of EigenSolver class and its subclasses
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                 \
        FUNC(test_PowerIteration) \
        FUNC(test_SlepcSolver)

#include "utilities/TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "callow/test/matrix_fixture.hh"
#include <iostream>

using namespace callow;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// system size
int n = 50;

// convergence
double tol    = 1e-10;
int    maxit  = 1000;

int test_PowerIteration(int argc, char *argv[])
{
  Vector X(n, 0.0);
  Vector X0(n, 1.0);
  Matrix::SP_matrix A = test_matrix_1(n);
  EigenSolver::SP_db db(new detran_utilities::InputDB("test_PowerIteration"));
  db->put<std::string>("eigen_solver_type", "power");
  db->put<double>("eigen_solver_tol", 1e-10);
  db->put<int>("eigen_solver_maxit", 1000);
  EigenSolver::SP_solver solver = EigenSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  solver->set_monitor_level(1);
  int status = solver->solve(X, X0);
  cout << solver->eigenvalue() << endl;
  return 0;
}

int test_SlepcSolver(int argc, char *argv[])
{
  Vector X(n, 0.0);
  Vector X0(n, 1.0);
  Matrix::SP_matrix A = test_matrix_1(n);
  EigenSolver::SP_db db(new detran_utilities::InputDB("test_SlepcSolver"));
  db->put<std::string>("eigen_solver_type", "slepc");
  db->put<double>("eigen_solver_tol", 1e-10);
  db->put<int>("eigen_solver_maxit", 1000);
  EigenSolver::SP_solver solver = EigenSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  solver->set_monitor_level(1);
  int status = solver->solve(X, X0);
  cout << solver->eigenvalue() << endl;
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
