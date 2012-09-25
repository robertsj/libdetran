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
// solvers
#include "callow/solver/PowerIteration.hh"
#ifdef CALLOW_ENABLE_SLEPC
#include "callow/solver/SlepcSolver.hh"
#endif
//
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
double abstol = 1e-13;
double reltol = 1e-13;
int    maxit  = 1000;

// reference
double X_ref[] = {5.459135698786395e+00, 8.236037378063020e+00,
                  9.648531187020394e+00, 1.036694341855932e+01,
                  1.073221653274062e+01, 1.091771229837346e+01,
                  1.101148972685055e+01, 1.105810668452673e+01,
                  1.107978760425382e+01, 1.108701173915058e+01,
                  1.108356356535519e+01, 1.106847349927400e+01,
                  1.103582874555248e+01, 1.097247463295948e+01,
                  1.085272174937750e+01, 1.062793308617560e+01,
                  1.020677234863001e+01, 9.418093128951856e+00,
                  7.941390927380783e+00, 5.176556370952343e+00};

int test_PowerIteration(int argc, char *argv[])
{
  Vector<double> X(n, 0.0);
  Vector<double> X0(n, 1.0);
  //X0[0] = 10;
  Matrix<double>::SP_matrix A = test_matrix_1<double>(n);
  A->print_matlab("A.out");
  PowerIteration<double> solver(1e-6, 10000);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_level(2);
  int status = solver.solve(X, X0);
//  TEST(status == 0);
//  for (int i = 0; i < 20; ++i)
//  {
//    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
//  }
  return 0;
}

int test_SlepcSolver(int argc, char *argv[])
{
#ifdef CALLOW_ENABLE_SLEPC
  Vector<double> X(n, 0.0);
  Vector<double> X0(n, 1.0);
  Matrix<double>::SP_matrix A = test_matrix_1<double>(n);
  A->print_matlab("A.out");
  SlepcSolver solver(1e-6, 10000);
  solver.set_operators(test_matrix_1<double>(n));
  solver.set_monitor_level(2);
  int status = solver.solve(X, X0);
  cout << solver.eigenvalue() << endl;
#endif
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
