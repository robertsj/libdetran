//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Track.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Track class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Track)

// Detran headers
#include "TestDriver.hh"
#include "Track.hh"

// Setup
/* ... */

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Track()
{
  // Test track creation.
  Point enter(0.0, 0.5);
  Point exit(0.5, 1.0);
  Track track(enter, exit);
  TEST(soft_equiv(track.enter().x(), 0.0));
  TEST(soft_equiv(track.exit().y(),  1.0));

  // Test track addition of segments.
  track.add_segment(Segment(0, 1.0));
  track.add_segment(Segment(1, 2.0));
  TEST(track.number_segments()   == 2);
  TEST(track.segment(0).region() == 0);
  TEST(track.segment(1).region() == 1);
  TEST(soft_equiv(track.segment(0).length(), 1.0));
  TEST(soft_equiv(track.segment(1).length(), 2.0));

  // Test iterators.
  Track::iterator it = track.begin();
  TEST((*it).region() == 0);
  it++;
  TEST((*it).region() == 1);
  Track::riterator rit = track.rbegin();
  TEST((*rit).region() == 1);
  rit++;
  TEST((*rit).region() == 0);

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Track.cc
//---------------------------------------------------------------------------//
