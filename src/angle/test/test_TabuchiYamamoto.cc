//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_TabuchiYamamoto.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of TabuchiYamamoto class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "TabuchiYamamoto.hh"

using namespace detran_utilities;
using namespace detran_angle;
using namespace std;

TEST(TabuchiYamamoto, Basic)
{
  // Construct quadrature.
  TabuchiYamamoto q(1);
  EXPECT_EQ(q.number_polar(), 1);
  EXPECT_NEAR(q.sin_theta(0),  0.798184, 1.0e-12);
  EXPECT_NEAR(q.weight(0),     1.000000, 1.0e-12);
}

//---------------------------------------------------------------------------//
//              end of test_TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
