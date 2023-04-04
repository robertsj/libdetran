//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Collocated.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Collocated class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Collocated.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(Collocated, Basic)
{
  // Construct quadrature.
  Collocated q(2, 5, 1, 2, "TY");
  EXPECT_EQ(q.number_angles()         , 4 * 5 * 2);
  EXPECT_EQ(q.number_angles_octant()  , 5 * 2);
  EXPECT_EQ(q.number_azimuths_octant(), 5);
  EXPECT_EQ(q.number_polar_octant()   , 2);
  EXPECT_EQ(q.number_tracks(0)        , 10);
  EXPECT_EQ(q.number_enter(0, 0)      , 1);
  EXPECT_EQ(q.number_enter(0, 1)      , 9);
  EXPECT_NEAR(q.enter(0, 0).x(), 0.0, 1.0e-12);
  EXPECT_NEAR(q.enter(0, 0).y(), 9.444444444444444e-01, 1.0e-12);

  //q.display();
  //q.display_tracks();
}

//---------------------------------------------------------------------------//
//              end of test_Collocated.cc                                    //
//---------------------------------------------------------------------------//
