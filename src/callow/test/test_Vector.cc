//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_Vector.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of Vector class.
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "callow/vector/Vector.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/Definitions.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// Test of basic public interface
TEST(Vector, Basic)
{
  {
    Vector v(10, 1.0);
    for (int i = 0; i < 10; i++)
    {
      EXPECT_NEAR(v[i], 1.0, 1.0e-12);
    }
    v[5] = 5.0;
    EXPECT_NEAR(v[5], 5.0, 1.0e-12);
    v[5] = 1.0;

    Vector y(10, 2.0);
    double val = v.dot(y);
    double val2 = 0.0;
    for (int i = 0; i < 10; i++)
    {
      val2 += v[i]*y[i];
    }
    cout << val2 << endl;
    EXPECT_NEAR(val, 20.0, 1.0e-12);

    v.set(1.23);
    EXPECT_NEAR(v[0], 1.23, 1.0e-12);
    v.scale(2.0);
    EXPECT_NEAR(v[0], 2.46, 1.0e-12);
  }

  {
    detran_utilities::vec_dbl v(5, 1.0);
    Vector V(v);
    EXPECT_EQ(V.size(), 5);
    EXPECT_NEAR(V[4], 1.0, 1.0e-12);
  }

  {
    double v[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    Vector V(5, v);
    EXPECT_EQ(V.size(), 5);
    EXPECT_NEAR(V[4], 1.0, 1.0e-12);
    double v_L1 = V.norm(L1);
    EXPECT_NEAR(v_L1, 5.0, 1.0e-12);
    double v_L2 = V.norm(L2);
    EXPECT_NEAR(v_L2, 2.236067977499790, 1.0e-12);
    double v_LI = V.norm(LINF);
    EXPECT_NEAR(v_LI, 1.0, 1.0e-12);
  } 
}

TEST(Vector, Resize)
{
  typedef Vector Vec_T;
  Vec_T v;
  v.resize(5, 1.0);
  v.set(3.0);
  v.display();
}

//---------------------------------------------------------------------------//
//              end of test_Vector.cc
//---------------------------------------------------------------------------//
