//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_NonlinearArnoldi.cc
 *  @brief Test of test_NonlinearArnoldi class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_NonlinearArnoldi)

#include "utilities/TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "callow/solver/NonlinearArnoldi.hh"
#include "callow/solver/PowerIteration.hh"

#include "matrix_fixture.hh"
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

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

class DiagonalPC: public PCShell
{
public:
  DiagonalPC(Matrix::SP_matrix A, Matrix::SP_matrix B, void* context)
    : PCShell("diag", context), d_A(A), d_B(B)
  {/* ... */}
  void apply(Vector &x, Vector &y)
  {
    std::cout << " applied" << std::endl;
    double e = ((NonlinearArnoldi*)PCShell::d_context)->eigenvalue();
    for (size_t i = 0; i < x.size(); ++i)
    {
      double val = (*d_A)[d_A->diagonal(i)] - 1 * (*d_B)[d_B->diagonal(i)];
      y[i] = x[i] / val;
    }
  }
private:
  Matrix::SP_matrix d_A;
  Matrix::SP_matrix d_B;
};

/*
 *  Solves Ax = e*B*x where A is the standard 2nd difference
 *  operator and B is diagonal matrix of 1 through n.
 */
int test_NonlinearArnoldi(int argc, char *argv[])
{
  Matrix::SP_matrix A = test_matrix_2(10);

  int n = A->number_columns();
  Vector X (n, 1.0);
  Vector X0(n, 1.0);

  A->print_matlab("A.out");
  // let's see if we can do Ax - eIx = 0 with it.  This means B = I.
  if (1) {
    double ref[] =
    { 1.000000000000000, -1.872738525887991, 2.554768633964161 };
    Matrix::SP_matrix B(new Matrix(n, n, 1));
    for (int i = 0; i < n; ++i) B->insert(i, i, 1);
    B->assemble();
    NonlinearArnoldi solver(1e-8, 50, 20);
    solver.set_operators(A, B);
    solver.solve(X, X0);
    X.scale(1.0/X[0]);
    for (int i = 0; i < 3; ++i)
    {
      printf("%16.8f  %16.8f \n", X[i], ref[i]);
      TEST(soft_equiv(X[i], ref[i], 1.0e-6));
    }
  }

  // let's see if we can do Ax - eIx = 0 with it.  This means B = I.
  if (0) {
    double ref[] =
    { 1.000000000000000, -1.872738525887991, 2.554768633964161 };
    Matrix::SP_matrix B(new Matrix(n, n, 2));
    for (int i = 0; i < n; ++i)
    {
      B->insert(i, i, 1.0);
      if (i < n-1) B->insert(i, i+1, 0.1);
    }
    B->assemble();
    B->display();
    NonlinearArnoldi solver(1e-8, 50, 50);
    PCShell::SP_preconditioner pc(new DiagonalPC(A, B, &solver));
    solver.set_operators(A, B);
    solver.set_preconditioner(pc);
    solver.solve(X, X0);
    X.scale(1.0/X[0]);
//    for (int i = 0; i < 3; ++i)
//      TEST(soft_equiv(X[i], ref[i], 1.0e-6));
  }

//  Matrix::SP_matrix B = new Matrix(n, n, 1);
//  for (int i = 0; i < n; ++i) B->insert(i, i, i+1);
//  B->assemble();
//  EigenSolver::SP_db db(new detran_utilities::InputDB("test_PowerIteration"));
//  db->put<std::string>("eigen_solver_type", "power");
//  db->put<double>("eigen_solver_tol", 1e-14);
//  db->put<int>("eigen_solver_maxit", 10000);
//  EigenSolver::SP_solver solver = new NonlinearArnoldi(1e-8, 50, 20);
//  solver->set_operators(test_matrix_1(n));
//  solver->set_monitor_level(2);
//  cout << solver->eigenvalue() << endl;
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
