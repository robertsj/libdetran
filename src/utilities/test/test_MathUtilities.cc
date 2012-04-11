/*
 * test_MathUtilities.cc
 *
 *  Created on: Apr 10, 2012
 *      Author: robertsj
 */


// LIST OF TEST FUNCTIONS
#define TEST_LIST            \
        FUNC(test_norm_L2)   \
        FUNC(test_norm_L1)   \
        FUNC(test_norm_Linf) \
        FUNC(test_vec_scale) \
        FUNC(test_norm_residual_L2) \
        FUNC(test_norm_residual_L1) \
        FUNC(test_norm_residual_Linf)


// Detran headers
#include "TestDriver.hh"
#include "MathUtilities.hh"

// System
#include <iostream>

using namespace detran_test;
using namespace detran;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

// Test definitions.

int test_norm_L2()
{
  vec_dbl v(10, 1.1);
  double norm_v = norm(v, "L2");
  TEST(soft_equiv(norm_v, 3.478505426185218));
  return 0;
}
int test_norm_L1()
{
  vec_dbl v(10, 1.1);
  double norm_v = norm(v, "L1");
  TEST(soft_equiv(norm_v, 11.0));
  return 0;
}
int test_norm_Linf()
{
  vec_dbl v(10, 1.1);
  double norm_v = norm(v, "Linf");
  TEST(soft_equiv(norm_v, 1.1));
  return 0;
}
int test_vec_scale()
{
  vec_dbl v(10, 1.1);
  vec_scale(v, 2.34);
  for (int i = 0; i < 10; i++)
  {
    TEST(soft_equiv(v[i], 2.574));
  }
  return 0;
}
int test_norm_residual_L2()
{
  vec_dbl v(10, 1.1);
  vec_dbl y(10, 1.0);
  double norm_res = norm_residual(v, y, "L2");
  TEST(soft_equiv(norm_res, 0.316227766016838));
  return 0;
}
int test_norm_residual_L1()
{
  vec_dbl v(10, 1.1);
  vec_dbl y(10, 1.0);
  double norm_res = norm_residual(v, y, "L1");
  TEST(soft_equiv(norm_res, 1.0));
  return 0;
}
int test_norm_residual_Linf()
{
  vec_dbl v(10, 1.1);
  vec_dbl y(10, 1.0);
  double norm_res = norm_residual(v, y, "Linf");
  TEST(soft_equiv(norm_res, 0.1));
  return 0;
}


//---------------------------------------------------------------------------//
//              end of testMathUtilities.cc
//---------------------------------------------------------------------------//




