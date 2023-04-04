//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Uniform.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Uniform class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Uniform.hh"

using namespace detran_utilities;
using namespace detran_angle;
using namespace std;

TEST(Uniform, Basic)
{
  // Construct quadrature.
  Uniform q(2, 5, 13, 3, "TY");
  EXPECT_EQ(q.number_angles(), 4 * 5 * 3);
  EXPECT_EQ(q.number_angles_octant(), 5 * 3);
  EXPECT_EQ(q.number_azimuths_octant(), 5);
  EXPECT_EQ(q.number_polar_octant(), 3);
  EXPECT_EQ(q.number_tracks(0), 13);
  EXPECT_EQ(q.number_enter(0, 0), 2);
  EXPECT_EQ(q.number_enter(0, 1), 11);
  EXPECT_EQ(q.number_exit(0, 0) , 11);
  EXPECT_EQ(q.number_exit(0, 1) , 2);
  double y = 1.0 - 0.5 / 11.0;
  EXPECT_NEAR(q.enter(0, 0).x(), 0.0, 1.0e-12);
  EXPECT_NEAR(q.enter(0, 0).y(), y, 1.0e-12);

  EXPECT_EQ(q.polar(13), 1);
  EXPECT_EQ(q.azimuth(13), 4);
  q.display_tracks();

  //
  Uniform q2(2, 1, 3, 1, "TY");
  EXPECT_EQ(q2.number_enter(0, 0), 1);
  EXPECT_EQ(q2.number_enter(0, 1), 2);
  EXPECT_EQ(q2.number_exit(0, 0) , 2);
  EXPECT_EQ(q2.number_exit(0, 1) , 1);
  q2.display_tracks();
}

//---------------------------------------------------------------------------//
//              end of test_TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
