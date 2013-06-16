//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Tracker.cc
 *  @brief Test of Tracker class
 *  @note  Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_Tracker_2x2)  \
        FUNC(test_Tracker_3x3)  \
        FUNC(test_Tracker_pin)

#include "TestDriver.hh"
#include "Tracker.hh"
#include "Mesh2D.hh"
#include "csg_fixture.hh"
#include "angle/QuadratureFactory.hh"
#include "ioutils/PSPlotter.hh"

using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Tracker_2x2(int argc, char *argv[])
{
//
//  // Create mesh
//  vec_dbl cm(2, 0.0);
//  cm[1] = 1.0;
//  vec_int fm(1, 2);
//  vec_int mat(1, 0);
//  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mat));
//
//  // Create quadrature
//  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));
//  quad->display_tracks();
//
//  // Create tracker
//  Tracker tracker(mesh, quad);
//  Tracker::SP_trackdb tracks = tracker.trackdb();
//
//  // Verify tracker.
//  double length = 0.559016994374947;
//  // Number of segments for each track within angle
//  int ns[] = {1,2,1};
//  // Region map for all segments
//  int region[] = {2,0,3,1,3,1,2,0};
//  // Region counter
//  int r = 0;
//  for (int a = 0; a < 2; a++)
//  {
//    TEST(tracks->number_tracks_angle(a) == 3);
//    for (int t = 0; t < 3; t++)
//    {
//      Tracker::SP_track track = tracks->track(a, t);
//      TEST(track->number_segments() == ns[t]);
//      for (int s = 0; s < ns[t]; s++)
//      {
//        TEST(soft_equiv(track->segment(s).length(), length));
//        TEST(track->segment(s).region() == region[r++]);
//      }
//    }
//  }
  return 0;
}

int test_Tracker_3x3(int argc, char *argv[])
{
//
//  // Create mesh
//  vec_dbl cm(2, 0.0);
//  cm[1] = 1.0;
//  vec_int fm(1, 3);
//  vec_int mat(1, 0);
//  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mat));
//  mesh->display();
//  // Create quadrature 2,1,3,1
//  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));
//  quad->display_tracks();
//
//  // Create tracker
//  Tracker tracker(mesh, quad);
//  Tracker::SP_trackdb tracks = tracker.trackdb();
//  tracks->display();
//
//  // Number of segments for each track within angle
//  int ns[] = {2,5,2};
//  // Region map for all segments
//  int region[] = {6,7, 0,3,4,5,8, 1,2,
//                  8,7, 2,5,4,3,6, 1,0};
//  double len0 = 0.186338998124982;
//  double len1 = 0.372677996249965;
//  double length[] = {len1,len0,  len0,len0,len1,len0,len0,  len0,len1};
//
//  // Region counter
//  int r = 0;
//  for (int a = 0; a < 2; a++)
//  {
//    TEST(tracks->number_tracks_angle(a) == 3);
//    // Length counter
//    int l = 0;
//    for (int t = 0; t < 3; t++)
//    {
//      Tracker::SP_track track = tracks->track(a, t);
//      TEST(track->number_segments() == ns[t]);
//      for (int s = 0; s < ns[t]; s++)
//      {
//        TEST(soft_equiv(track->segment(s).length(), length[l++]));
//        TEST(track->segment(s).region() == region[r++]);
//      }
//    }
//  }
//
//  // Normalize the lengths.
//  tracker.normalize();
//  tracks->display();
//  TEST(soft_equiv(tracks->track(0, 0)->segment(0).length(), 0.662538659999938));
  return 0;
}

// Test the tracking of a pin
int test_Tracker_pin(int argc, char *argv[])
{
  Geometry::SP_geometry pin = test_2D_pincell_simple();

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 0.3);
  db->put<std::string>("tracker_spatial_quad_type", "gl");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 3);

  Tracker::SP_quadrature q = detran_angle::QuadratureFactory::build(db, 2);

  Tracker tracker(db, q);

  tracker.trackit(pin);

  vec_dbl bbox(4, 0.0);
  bbox[1] = 1.27;
  bbox[3] = 1.27;
  detran_ioutils::PSPlotter plt("tracks.eps", bbox);
  Tracker::SP_track track = tracker.trackdb()->track(0, 0);
  plt.draw_line(track->enter(), track->exit());
  plt.finalize();

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Tracker.cc
//----------------------------------------------------------------------------//
