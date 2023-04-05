//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_QuadrupleRange.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of QuadrupleRange class
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "TestDriver.hh"
#include "AbuShumaysQuadrupleRange.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(QuadrupleRange, Basic)
{
  // Get quadrature fixture
  QuadrupleRange::SP_quadrature q = std::make_shared<QuadrupleRange>(2, 2, 1);
  EXPECT_TRUE(q != NULL);
  q->display();
  EXPECT_EQ(q->number_angles()  , 8);
  EXPECT_EQ(q->number_octants() , 4);
  EXPECT_NEAR(q->mu(0, 1),  0.2582870761957, 1.0e-12);
  EXPECT_NEAR(q->mu(0, 0),  0.7417129238043, 1.0e-12);
  EXPECT_NEAR(q->eta(0, 1), 0.7417129238043, 1.0e-12);
  EXPECT_NEAR(q->eta(0, 0), 0.2582870761957, 1.0e-12);
  EXPECT_NEAR(q->weight(0), 1.5707963267949, 1.0e-12);
  EXPECT_NEAR(q->weight(1), 1.5707963267949, 1.0e-12);
}

//---------------------------------------------------------------------------//
//              end of test_QuadrupleRange.cc
//---------------------------------------------------------------------------//
