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

#define COUT(c) std::cout << c << std::endl;


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
<<<<<<< HEAD
  InputDB::SP_input db =std::make_shared<InputDB>();
=======
  // initialize simple set of tracks
  InputDB::SP_input db = InputDB::Create();
>>>>>>> dev
  db->put<string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 2);
  db->put<int>("quad_number_polar_octant",   1);
  auto q = std::dynamic_pointer_cast<ProductQuadrature>(QF::build(db, 2));
  q->display();
  TrackDB::map_track all_tracks;
  for (int o = 0; o < 4; ++o)
  {
    for (int a = 0; a < 2; ++a)
    {
      TrackDB::vec_track tracks;
      for (int t = 0; t < a + 1; ++t)
      {
        Point entry(a, 0+t);
        Point exit(o, 1+t);
        auto track = std::make_shared<Track>(entry, exit, 1.0);
        Segment s(a, 1.0);
        track->add_segment(s);
        tracks.push_back(track);
      }
      all_tracks[{o, a}] = tracks;
    }
  }
  TrackDB trackdb(q, all_tracks);
  trackdb.display();

  // test sizes are correct
  TEST(trackdb.tracks().size() == 8);
  for (int o = 0; o < 4; ++o)
    for (int a = 0; a < 2; ++a)
      TEST(trackdb.tracks(o, a).size() == a+1);

  // test that the volume is correct
  auto volume = trackdb.volume();
  TEST(volume.size() == 2);
  double volume_0 = 0.5; // one track/segment of width 1 * length 1 * 0.5 (for half the angles)
  double volume_1 = 1.0; // two track/segments of width 1 * length 1 * 0.5 (for half the angles)

  COUT("vol 0 = " << volume[0]);
  COUT("vol 1 = " << volume[1]);

<<<<<<< HEAD
//----------------------------------------------------------------------------//
int test_TrackDB_3D(int argc, char *argv[])
{
  InputDB::SP_input db =std::make_shared<InputDB>();
  db->put<string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);
  db->put<int>("quad_number_polar_octant", 2);
=======
  TEST(volume[0] == volume_0);
  TEST(volume[1] == volume_1);
>>>>>>> dev


  // test normalization.  now a=0 segments should have length 4.2
  //  and a=1 segments should have length 2.1
  std::vector<double> ref_vols(2, 2.1);
  trackdb.normalize(ref_vols);

  // first check the volume is updated
  volume = trackdb.volume();

  COUT("vol 0 = " << volume[0]);
  COUT("vol 1 = " << volume[1]);

  TEST(soft_equiv(volume[0], 2.1));
  TEST(soft_equiv(volume[1], 2.1));

  TEST(soft_equiv(trackdb.tracks(0, 0)[0]->segment(0).length(), 4.2));
  TEST(soft_equiv(trackdb.tracks(0, 1)[0]->segment(0).length(), 2.1));

  trackdb.display();

  return 0;
}

//----------------------------------------------------------------------------//
int test_TrackDB_3D(int argc, char *argv[])
{
  // Nothing here till 3D tracking is implemented.

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TrackDB.cc
//---------------------------------------------------------------------------//
