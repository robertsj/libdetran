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
  Vector<double> X(n, 0.0);
  Vector<double> X0(n, 1.0);
  Matrix<double>::SP_matrix A = test_matrix_1<double>(n);
  EigenSolver<double>::SP_solver solver =
    EigenSolverCreator<double>::Create("power", tol, maxit);
  solver->set_operators(test_matrix_1<double>(n));
  solver->set_monitor_level(1);
  int status = solver->solve(X, X0);
  cout << solver->eigenvalue() << endl;
  return 0;
}

int test_SlepcSolver(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> X0(n, 1.0);
  Matrix<double>::SP_matrix A = test_matrix_1<double>(n);
  EigenSolver<double>::SP_solver solver =
    EigenSolverCreator<double>::Create("slepc", tol, maxit);
  solver->set_operators(test_matrix_1<double>(n));
  solver->set_monitor_level(1);
  int status = solver->solve(X, X0);
  cout << solver->eigenvalue() << endl;
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
