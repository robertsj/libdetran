//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_LinearSolver.cc
 *  @brief Test of LinearSolver class and its subclasses
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST              \
        FUNC(test_Richardson)  \
        FUNC(test_Jacobi)      \
        FUNC(test_GaussSeidel) \
        FUNC(test_SOR)         \
        FUNC(test_GMRES)       \
        FUNC(test_PetscSolver)

#include "utilities/TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "LinearSolverCreator.hh"
// solvers
#include "callow/solver/Richardson.hh"
#include "callow/solver/Jacobi.hh"
#include "callow/solver/GaussSeidel.hh"
#include "callow/solver/GMRES.hh"
#ifdef CALLOW_ENABLE_PETSC
#include "callow/solver/PetscSolver.hh"
#endif
// pc
#include "callow/preconditioner/PCJacobi.hh"
#include "callow/preconditioner/PCILU0.hh"
#include "callow/preconditioner/PCIdentity.hh"
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
int n = 20;
// solver
LinearSolverCreator::SP_solver solver;
// database
LinearSolverCreator::SP_db db;

LinearSolver::SP_db get_db()
{
  LinearSolver::SP_db p(new detran_utilities::InputDB("callow_db"));
  p->put<int>("linear_solver_maxit",   1000);
  p->put<double>("linear_solver_atol", 1e-13);
  p->put<double>("linear_solver_rtol", 1e-13);
  p->put<int>("linear_solver_monitor_level", 2);
  p->put<int>("linear_solver_monitor_diverge", 0);
  return p;
}

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

int test_Richardson(int argc, char *argv[])
{
  Vector X(n, 0.0);
  Vector B(n, 1.0);
  db = get_db();
  db->put<std::string>("linear_solver_type", "richardson");
  //db->put<int>("linear_solver_monitor_level", 1);
  solver = LinearSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  //Preconditioner::SP_preconditioner pcilu0(new PCILU0(test_matrix_1(n)));
  //solver->set_preconditioner(pcilu0, LinearSolver::LEFT);
  int status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_Jacobi(int argc, char *argv[])
{
  Vector X(n, 0.0);
  Vector B(n, 1.0);
  db = get_db();
  db->put<std::string>("linear_solver_type", "jacobi");
  solver = LinearSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  int status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_GaussSeidel(int argc, char *argv[])
{
  Vector X(n, 0.0);
  Vector B(n, 1.0);
  db = get_db();
  db->put<std::string>("linear_solver_type", "gauss-seidel");
  solver = LinearSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  int status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_SOR(int argc, char *argv[])
{
  Vector X(n, 0.0);
  Vector B(n, 1.0);
  db = get_db();
  db->put<std::string>("linear_solver_type", "sor");
  db->put<double>("linear_solver_sor_omega", 1.3);
  solver = LinearSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  int status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
  return 0;
}

int test_GMRES(int argc, char *argv[])
{

  GMRES::SP_matrix A;
  A = test_matrix_1(n);
  Vector X(n, 0.0);
  Vector B(n, 1.0);
  db = get_db();
  db->put<std::string>("linear_solver_type", "gmres");
  db->put<int>("linear_solver_maxit", 50);
  db->put<int>("linear_solver_gmres_restart", 16);
  // NO PC
  std::cout << "*** GMRES + NO PC ***" << std::endl;
  solver = LinearSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n));
  int status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }

  Preconditioner::SP_preconditioner pcilu0;
  Preconditioner::SP_preconditioner pcjacobi;
  pcilu0 = new PCILU0(A);
  pcjacobi = new PCJacobi(A);

  // PCILU0 -- LEFT
  std::cout << "*** GMRES + ILU(0) on LEFT ***" << std::endl;
  X.set(0.0);
  solver->set_preconditioner(pcilu0, GMRES::LEFT);
  status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }

  // PCILU0 -- RIGHT
  std::cout << "*** GMRES + ILU(0) on RIGHT ***" << std::endl;
  X.set(0.0);
  solver->set_preconditioner(pcilu0, GMRES::RIGHT);
  status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }

  // PCJacobi -- LEFT
  std::cout << "*** GMRES + Jacobi on LEFT ***" << std::endl;
  X.set(0.0);
  solver->set_preconditioner(pcjacobi, GMRES::LEFT);
  status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }

  // PCJacobi -- RIGHT
  std::cout << "*** GMRES + Jacobi on RIGHT ***" << std::endl;
  X.set(0.0);
  solver->set_preconditioner(pcjacobi, GMRES::RIGHT);
  status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }

  return 0;
}

int test_PetscSolver(int argc, char *argv[])
{
#ifdef CALLOW_ENABLE_PETSC
  Vector X(n, 0.0);
  Vector B(n, 1.0);
  db = get_db();
  db->put<std::string>("linear_solver_type", "petsc");
  db->put<std::string>("pc_type", "ilu0");
//  db->put<std::string>("pc_type", "petsc_pc");
//  db->put<std::string>("petsc_pc_type", "ilu");
  solver = LinearSolverCreator::Create(db);
  solver->set_operators(test_matrix_1(n), db);
  int status = solver->solve(B, X);
  TEST(status == SUCCESS);
  for (int i = 0; i < 20; ++i)
  {
    TEST(soft_equiv(X[i],  X_ref[i], 1e-9));
  }
#endif
  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_LinearSolver.cc
//---------------------------------------------------------------------------//
