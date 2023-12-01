//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Threading.cc
 *  @brief Test of OpenMP threading in callow
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "callow/utils/Initialization.hh"
// threaded objects
#include "callow/solver/Jacobi.hh"
#include "callow/matrix/Matrix.hh"
//
#include "callow/test/matrix_fixture.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

#define COUT(c) cout << c << endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
TEST(Threading, ThreadedMatVec)
{
    {
    int n = 200;
    Matrix::SP_matrix A = test_matrix_2(n);
    Vector x(A->number_columns(), 1.0);
    Vector y(A->number_columns(), 0.0);
    int nt = 30; // number of trials
    int i = 0;
    //double t = TIME;
    #pragma omp parallel default(shared) private(i)
    {
      for (i = 0; i < nt; ++i)
      {
        A->multiply(x, y);
      }
    }
    COUT("y[0:3]=" << y[0] << " " << y[1] << " " << y[2])
    COUT("norm y = " << y.norm())
   // t = omp_get_wtime() - t;
    //COUT("ELAPSED = " << t)
  }
  COUT("done!!")
  callow_finalize();
}

//----------------------------------------------------------------------------//
TEST(Threading, ThreadedJacobi)
{
    {
    int n = 100;
    Matrix::SP_matrix A = test_matrix_2(n);
    Vector x(A->number_columns(), 0.0);
    Vector b(A->number_columns(), 1.0);
   // double t = omp_get_wtime();
    Jacobi solver(1.0e-8, 1.0e-8, 10000, 1.0);
    solver.set_operators(A);
    solver.set_monitor_level(1);
    #pragma omp parallel default(shared)
    {
      int status = solver.solve(b, x);
    }
    //t = omp_get_wtime() - t;
    //COUT("ELAPSED = " << t)
  }
  callow_finalize();
}

//----------------------------------------------------------------------------//
//              end of test_Threading.cc
//----------------------------------------------------------------------------//
