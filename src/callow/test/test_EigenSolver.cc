//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_EigenSolver.cc
 *  @brief Test of EigenSolver class and its subclasses
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "callow/utils/Initialization.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "callow/solver/SlepcSolver.hh"
#include "callow/preconditioner/PCJacobi.hh"
#include "callow/preconditioner/PCILU0.hh"
#include "callow/preconditioner/PCIdentity.hh"
#include "callow/test/matrix_fixture.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// system size
int n = 50;

// convergence
double tol    = 1e-10;
int    maxit  = 1000;

//----------------------------------------------------------------------------//
TEST(EigenSolver, PowerIteration)
{
  Vector X(n, 0.0);
  Vector X0(n, 1.0);
  Matrix::SP_matrix A = test_matrix_1(n);
  EigenSolver::SP_db db(new detran_utilities::InputDB("test_PowerIteration"));
  db->put<std::string>("eigen_solver_type", "power");
  db->put<double>("eigen_solver_tol", 1e-14);
  db->put<int>("eigen_solver_maxit", 10000);
  EigenSolver::SP_solver solver = EigenSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  solver->set_monitor_level(2);
  cout << solver->eigenvalue() << endl;
}

//----------------------------------------------------------------------------//
TEST(EigenSolver, SlepcSolver)
{
  Vector X(n, 0.0);
  Vector X0(n, 1.0);
  Matrix::SP_matrix A = test_matrix_1(n);
  EigenSolver::SP_db db(new detran_utilities::InputDB("test_SlepcSolver"));
  db->put<std::string>("eigen_solver_type", "slepc");
  db->put<double>("eigen_solver_tol", 1e-10);
  db->put<int>("eigen_solver_maxit", 10000);
  EigenSolver::SP_solver solver = EigenSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  solver->set_monitor_level(1);
  cout << solver->eigenvalue() << endl;
}

#include "matrix/MatrixShell.hh"
class MatShellGD: public MatrixShell
{
public:
  MatShellGD(MatrixBase::SP_matrix A)
    : MatrixShell(this, A->number_columns(), A->number_columns())
    , d_B(A)
  {
    /* ... */
  }
  void multiply(const Vector &x,  Vector &y)
  {
    d_B->multiply(x, y);
  }
  void multiply_transpose(const Vector &x,  Vector &y)
  {
    d_B->multiply_transpose(x, y);
  }
private:
  MatrixBase::SP_matrix d_B;
};

//----------------------------------------------------------------------------//
TEST(EigenSolver, SlepcGD)
{
  // Solve a generalized eigenvalue problem with GD.
  //   A * x = e * B * x
  // -->  (A-eB)x = 0
  //
  // the preconditioner M should approximate  (A-eB), so its
  // application should be ~ inv(A-eB).  Here, we'll set M = A-eB but
  // apply it as an ILU pc.
  int m = 10;
  Matrix::SP_matrix A = test_matrix_2(m);
  m = A->number_columns();
  Matrix::SP_matrix B = std::make_shared<Matrix>(m, m, 1);
  for (int i = 0; i < m; ++i) B->insert(i, i, 0.176376589713020);
  B->assemble();
  A->print_matlab("A.out");
  B->print_matlab("B.out");

  Vector X1(m, 0.0);
  Vector X0(m, 1.0);

  Matrix::SP_matrix C(new Matrix(*A));
  int    *D = C->diagonals();
  double *V = C->values();
  for (int i = 0; i < m; ++i)
    V[D[i]] -= 0.176376589713020;
  C->print_matlab("C.out");
  EigenSolver::SP_preconditioner P(new PCILU0(C));

  MatrixBase::SP_matrix S(new MatShellGD(C));

  SlepcSolver::SP_db db(new detran_utilities::InputDB("test_SlepcSolver"));
  db->put<std::string>("eigen_solver_type", "slepc");
  db->put<std::string>("eigen_solver_slepc_type", "gd");
  db->put<double>("eigen_solver_tol", 1e-5);
  db->put<int>("eigen_solver_maxit",  100);

  SlepcSolver solver("gd", 1.0e-5, 100, 1);

  solver.set_operators(A, B);
  solver.set_preconditioner_matrix(S);
  solver.set_monitor_level(1);
  solver.solve(X1, X0);
  cout << " EV = " << solver.eigenvalue() << endl;
}

//----------------------------------------------------------------------------//
//              end of test_EigenSolver.cc
//----------------------------------------------------------------------------//
