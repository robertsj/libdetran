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

TEST(Tracker2, TrackerBasic)
{
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1,  2);
  vec_int mt(1,  0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 1.0);
  db->put<int>("tracker_minimum_number", 3);
  db->put<std::string>("tracker_spatial_quad_type", "uniform");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);

  auto q = std::dynamic_pointer_cast<ProductQuadrature>(QuadratureFactory::build(db, 2));

  Tracker tracker(db, q);
  auto tracks = tracker.trackit(mesh);
  //tracker.trackdb()->display();
  //for (auto track: tracker.tracks())
  //  COUT(*track);



  // Segment lengths
  double a = sqrt(2.0)/6.0;
  double b = sqrt(2.0)/3.0;
  double c = sqrt(2.0)/2.0;
  double lengths[] = {a, c, b,a,b, b,a,b, c, a,
                      a, c, b,a,b, b,a,b, c, a};
  // Number of segments for each track within angle
  int ns[] = {1,1,3,3,1,1,  1,1,3,3,1,1};
  // Region map for all segments
  int region[] = {2, 2, 0,2,3, 0,1,3, 1, 1,
                  0, 0, 1,0,2, 1,3,2, 3, 3};


}
