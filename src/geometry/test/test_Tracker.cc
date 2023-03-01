#include <gtest/gtest.h>

#include "Tracker.hh"
#include "RegionFactory.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "csg_fixture.hh"
#include "angle/QuadratureFactory.hh"
#include "angle/UserQuadrature.hh"
#include "ioutils/PSPlotter.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>


using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace std;
#define COUT(c) std::cout << c << std::endl;


// Checks tracking of 2-d mesh
TEST(Tracker, TrackerBasic2D)
{
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1,  2);
  vec_int mt(1,  0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));


  /*       1/6    1/2     5/7
   *     ----------|----------
   *     |  /    /
   *  5/6|/    /   |
   *     |   /     |
   *     | /     / |
   *  1/2|---------|
   *     |     /
   *     |    /
   *     |  /
   *  1/6|/
   *     |------------
   */

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 1.0);
  db->put<int>("tracker_minimum_number", 3);
  db->put<std::string>("tracker_spatial_quad_type", "uniform");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);

  auto q = std::dynamic_pointer_cast<ProductQuadrature>(QuadratureFactory::build(db, 2));

  Tracker tracker(db, q);
  auto trackdb = tracker.trackit(mesh);
  trackdb->display();

  // Segment lengths
  double a = sqrt(2.0)/6.0; // 0.2357
  double b = sqrt(2.0)/3.0; // 0.4714
  double c = sqrt(2.0)/2.0; // 0.7071

  auto tracks = trackdb->tracks(0, 0);
  EXPECT_EQ(tracks.size(), 6);
  EXPECT_NEAR(distance(tracks[0]->enter(), Point(0.0, 1.0/6.0)), 0.0, 1e-12);
  EXPECT_NEAR(distance(tracks[0]->exit(),  Point(5.0/6.0, 1.0)), 0.0, 1e-12);
  EXPECT_EQ(tracks[0]->number_segments(), 3);
  EXPECT_DOUBLE_EQ(tracks[0]->segment(0).length(), b);
  EXPECT_DOUBLE_EQ(tracks[0]->segment(1).length(), a);
  EXPECT_DOUBLE_EQ(tracks[0]->segment(2).length(), b);

  EXPECT_NEAR(distance(tracks[1]->enter(), Point(0.0, 0.5)), 0.0, 1e-12);
  EXPECT_NEAR(distance(tracks[1]->exit(),  Point(0.5, 1.0)), 0.0, 1e-12);
  EXPECT_EQ(tracks[1]->number_segments(), 1);
  EXPECT_DOUBLE_EQ(tracks[1]->segment(0).length(), c);

  // \todo Probably check all quadrants!
}
