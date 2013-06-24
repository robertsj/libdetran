//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  test_TrackDB.cc
 *  @brief Test of TrackDB class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST             \
        FUNC(test_TrackDB_2D) \
        FUNC(test_TrackDB_3D)

#include "TestDriver.hh"
#include "TrackDB.hh"
#include "angle/QuadratureFactory.hh"

using namespace detran_angle;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

typedef QuadratureFactory QF;

int test_TrackDB_2D(int argc, char *argv[])
{
  InputDB::SP_input db = InputDB::Create();
  db->put<string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 2);
  db->put<int>("quad_number_polar_octant",   1);
  QF::SP_quadrature q = QF::build(db, 2);

  TrackDB trackdb(q);

  // Add tracks for only the first two quadrants; we go in reverse otherwise.
  for (int o = 0; o < 2; ++o)
  {
    for (int a = 0; a < 2; ++a)
    {
      for (int t = 0; t < a + 1; ++t)
      {
      TrackDB::SP_track track(new Track(Point(a, 0+t), Point(1-a, 1+t), 1.0));
      trackdb.add_track(a + 2*o, 0, track);
      }
    }
  }
  // Test basics
  TEST(trackdb.number_tracks(0) == 1);
  TEST(trackdb.number_tracks(1) == 2);
  TEST(trackdb.number_tracks(2) == 1);
  TEST(trackdb.number_tracks(3) == 2);

  // Go around all four quadrants and iterate over the tracks in each one
  for (int o = 0; o < 4; ++o)
  {
    for (int a = 0; a < 2; ++a)
    {
      int az = a + 2*o;
      TrackDB::iterator_angle track = trackdb.begin(az, 0);
      for (; track != trackdb.end(az, 0); ++track)
      {
        std::cout << " azimuth= " << az
                  << " entrance=" << (*track)->enter() << std::endl;
      }
    }
  }

  // Iterate over all tracks at once
//  TrackDB::iterator track = trackdb.begin();
//  for (; track != trackdb.end(); ++track)
//  {
//    std::cout << " entrance=" << track->enter() << std::endl;
//  }

  return 0;
}

//----------------------------------------------------------------------------//
int test_TrackDB_3D(int argc, char *argv[])
{
  InputDB::SP_input db = InputDB::Create();
  db->put<string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);
  db->put<int>("quad_number_polar_octant", 2);
  QF::SP_quadrature q = QF::build(db, 3);

  TrackDB tracks(q);

  for (int a = 0; a < 4; ++a)
  {
    for (int p = 0; p < 2; ++p)
    {
      TrackDB::SP_track t(new Track(Point(a, 0), Point(a+1, 1+p), 1.0));
      tracks.add_track(a, p, t);
    }
  }

  TEST(tracks.number_tracks(0, 0) == 1);
  TEST(tracks.number_tracks(0, 1) == 1);
  TEST(tracks.number_tracks(1, 0) == 1);
  TEST(tracks.number_tracks(1, 1) == 1);
  TEST(tracks.number_tracks(2, 0) == 1);
  TEST(tracks.number_tracks(2, 1) == 1);
  TEST(tracks.number_tracks(3, 0) == 1);
  TEST(tracks.number_tracks(3, 1) == 1);


  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TrackDB.cc
//---------------------------------------------------------------------------//
