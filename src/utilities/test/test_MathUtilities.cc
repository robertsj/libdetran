/*
 * test_MathUtilities.cc
 *
 *  Created on: Apr 10, 2012
 *      Author: robertsj
 */

#include <gtest/gtest.h>

// LIST OF TEST FUNCTIONS
#define TEST_LIST            \
        FUNC(test_norm_L2)   \
        FUNC(test_norm_L1)   \
        FUNC(test_norm_Linf) \
        FUNC(test_vec_scale) \
        FUNC(test_norm_residual_L2) \
        FUNC(test_norm_residual_L1) \
        FUNC(test_norm_residual_Linf) \
        FUNC(test_linspace) \
        FUNC(test_range)


// Detran headers
#include "MathUtilities.hh"

// System
#include <iostream>

using namespace detran_utilities;
using namespace std;


// Test definitions.

TEST(MathUtilities, NormL2)
{
  vec_dbl v(10, 1.1);
  double norm_v = norm(v, "L2");
  EXPECT_NEAR(norm_v, 3.478505426185218, 1.0e-12);
}
TEST(MathUtilities, NormL1)
{
  vec_dbl v(10, 1.1);
  double norm_v = norm(v, "L1");
  EXPECT_NEAR(norm_v, 11.0, 1.0e-12);
}
TEST(MathUtilities, NormLinf)
{
  vec_dbl v(10, 1.1);
  double norm_v = norm(v, "Linf");
  EXPECT_NEAR(norm_v, 1.1, 1.0e-12);
}
TEST(MathUtilities, VecScale)
{
  vec_dbl v(10, 1.1);
  vec_scale(v, 2.34);
  for (int i = 0; i < 10; i++)
  {
    EXPECT_NEAR(v[i], 2.574, 1.0e-12);
  }
}
TEST(MathUtilities, NormResidualL2)
{
  vec_dbl v(10, 1.1);
  vec_dbl y(10, 1.0);
  double norm_res = norm_residual(v, y, "L2");
  EXPECT_NEAR(norm_res, 0.316227766016838, 1.0e-12);
}
TEST(MathUtilities, NormResidualL1)
{
  vec_dbl v(10, 1.1);
  vec_dbl y(10, 1.0);
  double norm_res = norm_residual(v, y, "L1");
  EXPECT_NEAR(norm_res, 1.0, 1.0e-12);
}
TEST(MathUtilities, NormResidualLinf)
{
  vec_dbl v(10, 1.1);
  vec_dbl y(10, 1.0);
  double norm_res = norm_residual(v, y, "Linf");
  EXPECT_NEAR(norm_res, 0.1, 1.0e-12);
}
TEST(MathUtilities, LinSpace)
{
  vec_dbl v = linspace(0, 4, 5);
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    std::cout << v[i] << " " << i << std::endl;
    EXPECT_NEAR(v[i], double(i), 1.0e-12);
  }
  
  v = linspace(-2, 2, 5);
  EXPECT_NEAR(v[0], -2.0, 1.0e-12);
  EXPECT_NEAR(v[1], -1.0, 1.0e-12);
  
  double ref[] = {0.4, 1.2, 2.0, 2.8, 3.6};
  v = linspace_center(0, 4, 5);
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    std::cout << v[i] << " " << ref[i] << std::endl;
    EXPECT_NEAR(v[i], ref[i], 1.0e-12);
  }
  
  v = linspace_center(-2, 2, 4); // -1.5 -0.5 0.5 1.5
  EXPECT_NEAR(v[0], -1.5, 1.0e-12);
  EXPECT_NEAR(v[1], -0.5, 1.0e-12);
}

TEST(MathUtilities, Range)
{
  vec_int v = range<int>(0, 4);
  EXPECT_EQ(v.size(), 4);
  EXPECT_EQ(v[0], 0);
  EXPECT_EQ(v[1], 1);
  EXPECT_EQ(v[2], 2);
  EXPECT_EQ(v[3], 3);

  v = range<int>(4, 0);
  EXPECT_EQ(v.size(), 4);
  EXPECT_EQ(v[0], 4);
  EXPECT_EQ(v[1], 3);
  EXPECT_EQ(v[2], 2);
  EXPECT_EQ(v[3], 1);

  v = range<int>(0, 4, true);
  EXPECT_EQ(v.size(), 5);
  EXPECT_EQ(v[0], 0);
  EXPECT_EQ(v[1], 1);
  EXPECT_EQ(v[2], 2);
  EXPECT_EQ(v[3], 3);
  EXPECT_EQ(v[4], 4);

  v = range<int>(0, 0);
  EXPECT_EQ(v.size(), 0);

}


//---------------------------------------------------------------------------//
//              end of test_MathUtilities.cc
//---------------------------------------------------------------------------//




