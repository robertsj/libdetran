//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Track.cc
 *  @brief Test of Track class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Track.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;

TEST(Track, Basic)
{
  // Test track creation.
  Point enter(0.0, 0.5);
  Point exit(0.5, 1.0);
  Track track(enter, exit, 0.123);
  EXPECT_NEAR(track.enter().x(), 0.0, 1.0e-12);
  EXPECT_NEAR(track.exit().y(),  1.0, 1.0e-12);
  EXPECT_NEAR(track.exit().z(),  0.0, 1.0e-12);
  EXPECT_NEAR(track.spatial_weight(),     0.123, 1.0e-12);

  // Test track addition of segments.
  track.add_segment(Segment(0, 1.0));
  track.add_segment(Segment(1, 2.0));
  EXPECT_EQ(track.number_segments()  , 2);
  EXPECT_EQ(track.segment(0).region(), 0);
  EXPECT_EQ(track.segment(1).region(), 1);
  EXPECT_NEAR(track.segment(0).length(), 1.0, 1.0e-12);
  EXPECT_NEAR(track.segment(1).length(), 2.0, 1.0e-12);

  // Test iterator
  std::cout << track << std::endl;
  Track::iterator it = track.begin();
  std::cout << " region = " << it->region() << std::endl;

  EXPECT_EQ((*it).region(), 0);
  it++;
  EXPECT_EQ((*it).region(), 1);

  // Test reverse
  Track::SP_track p_rtrack = track.reverse();
  auto &rtrack = *p_rtrack;
  std::cout << rtrack << std::endl;
  it = rtrack.begin();
  *it;
  EXPECT_EQ((*it).region(), 1);
  it++;
  EXPECT_EQ((*it).region(), 0);
}

//---------------------------------------------------------------------------//
//              end of test_Track.cc
//---------------------------------------------------------------------------//
