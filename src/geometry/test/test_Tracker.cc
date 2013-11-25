//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Tracker.cc
 *  @brief Test of Tracker class
 *  @note  Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Tracker_2d_mesh)    \
        FUNC(test_Tracker_3d_mesh)    \
        FUNC(test_Tracker_pin_2d)     \
        FUNC(test_Tracker_box_3d)

#include "TestDriver.hh"
#include "Tracker.hh"
#include "RegionFactory.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "csg_fixture.hh"
#include "angle/QuadratureFactory.hh"
#include "ioutils/PSPlotter.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;
#define COUT(c) std::cout << c << std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

typedef QuadraticSurfaceFactory   QSF;

int test_Tracker_2d_mesh(int argc, char *argv[])
{
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1,  2);
  vec_int mt(1,  0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 1.0);
  db->put<std::string>("tracker_spatial_quad_type", "uniform");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);

  Tracker::SP_quadrature q = detran_angle::QuadratureFactory::build(db, 2);

  Tracker tracker(db, q);
  tracker.trackit(mesh);
  //tracker.trackdb()->display();

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

  int i = 0;
  for (int a = 0; a < 2; ++a)
  {
    TEST(tracker.trackdb()->number_tracks(a, 0) == 6);
    for (int t = 0; t < tracker.trackdb()->number_tracks(a, 0); ++t)
    {
      Tracker::SP_track track = tracker.trackdb()->track(a, 0, t);
      //COUT(*track)
      TEST(track->number_segments() == ns[t]);
      for (int s = 0; s < track->number_segments(); ++s, ++i)
      {
        TEST(soft_equiv(track->segment(s).length(), lengths[i]));
        //COUT(track->segment(s).region() << " " << region[i])
        TEST(track->segment(s).region() == region[i]);
      }
    }
  }

  return 0;
}

//----------------------------------------------------------------------------//
int test_Tracker_3d_mesh(int argc, char *argv[])
{
  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1,  2);
  vec_int mt(1,  0);
  Mesh::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 1.0);
  db->put<std::string>("tracker_spatial_quad_type", "uniform");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);
  db->put<int>("quad_number_polar_octant",   1);

  Tracker::SP_quadrature q = detran_angle::QuadratureFactory::build(db, 3);
  q->display();
  // 0     0.6123724356958    0.6123724356958    0.5000000000000    1.5707963267949

  Tracker tracker(db, q);

  Tracker::SP_geometry geo = tracker.mesh_to_geometry(mesh);

  // enter=(0.000000, 0.833333, 0.166667) exit=(0.166667, 1.000000, 0.302749)
  double z = 1. / (12. * std::cos(0.25*std::acos(-1.))*sqrt(1.-0.25)) + 1./6.;

  Point P0(0.0, 5.0/6.0, 1.0/6.0);
  Point P1(1./6., 1., z);

  Point D(q->mu(0, 0), q->eta(0, 0), q->xi(0, 0));
  Ray R(P0, D);
  Tracker::SP_region r = geo->region(0);
  COUT(r->bound_max() << r->bound_min());
  bool v = geo->region(0)->intersects_bounding_box(R, 10.0);
  COUT("v fucking is " << v)
  tracker.trackit(mesh);
  tracker.trackdb()->display();

  return 0;

  double a = sqrt(3.0)/6.0;
  double b = sqrt(3.0)/3.0;
  double c = sqrt(3.0)/2.0;
  double lengths[] = {a, c, b, a, b, b, a, b, c, a,
                      a, c, b, a, b, b, a, b, c, a};

  // Number of segments for each track within angle
  int ns[] = {1,1,3,3,1,1,  1,1,3,3,1,1};
  // Region map for all segments
  int region[] = {2, 2, 0, 2, 3, 0, 1, 3, 1, 1,
                  0, 0, 1, 0, 2, 1, 3, 2, 3, 3};

  for (int a = 0; a < 4; ++a)
  {
    for (int p = 0; p < 2; ++a)
    {
      TEST(tracker.trackdb()->number_tracks(a, 0) == 6);
      int i = 0;
      for (int t = 0; t < tracker.trackdb()->number_tracks(a, p); ++t)
      {
        Tracker::SP_track track = tracker.trackdb()->track(a, p, t);
        TEST(track->number_segments() == ns[t]);
        for (int s = 0; s < track->number_segments(); ++s, ++i)
        {
          TEST(soft_equiv(track->segment(s).length(), lengths[i]));
          TEST(track->segment(s).region() == region[i]);
        }
      }
    }
  }

  return 0;
}

//----------------------------------------------------------------------------//
int test_Tracker_pin_2d(int argc, char *argv[])
{
  // P = 1.26, R = 0.54 --> V_f = 0.916088417786784, V_m =  0.0839115822132163

  Geometry::SP_geometry pin = test_2D_pincell_simple();

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 0.7);
  db->put<std::string>("tracker_spatial_quad_type", "uniform");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);

  Tracker::SP_quadrature q = detran_angle::QuadratureFactory::build(db, 2);

  COUT("FUCK THIS BOX: " << pin->region(1)->bound_min() << " " << pin->region(1)->bound_max())
//  Tracker tracker(db, q);
//
//  tracker.trackit(pin);
//
//  tracker.trackdb()->display();


  return 0;
}

// Track something shifted
int test_Tracker_box_3d(int argc, char *argv[])
{


  typedef Surface::SP_surface       SP_surface;
  typedef Surface::vec_surface      vec_surface;
  typedef Region::SP_node           SP_node;
  typedef std::vector<SP_node>      vec_node;
  typedef Region::SP_region         SP_region;
  typedef CSG_Node::vec_point       vec_point;

  vec_surface surfaces(6);
  surfaces[0] = QSF::CreatePlaneX(0.0);
  surfaces[1] = QSF::CreatePlaneX(1.0);
  surfaces[2] = QSF::CreatePlaneY(0.0);
  surfaces[3] = QSF::CreatePlaneY(1.0);
  surfaces[4] = QSF::CreatePlaneZ(0.0);
  surfaces[5] = QSF::CreatePlaneZ(1.0);

  SP_region box = Region::Create(0, Point(0,0,0)-0.001, Point(1,1,1)+.001);
  box->append(surfaces[0], true);
  box->append(surfaces[1], false);
  box->append(surfaces[2], true);
  box->append(surfaces[3], false);
  box->append(surfaces[4], true);
  box->append(surfaces[5], false);

  Geometry::SP_geometry geo = Geometry::Create(1.0, 1.0, 1.0);
  geo->add_region(box);

  InputDB::SP_input db = InputDB::Create();
  db->put<double>("tracker_maximum_spacing", 0.1);
  db->put<std::string>("tracker_spatial_quad_type", "uniform");
  db->put<std::string>("quad_type", "u-dgl");
  db->put<int>("quad_number_azimuth_octant", 1);
  db->put<int>("quad_number_polar_octant", 1);

  Tracker::SP_quadrature q = detran_angle::QuadratureFactory::build(db, 3);

  Tracker tracker(db, q);

  tracker.trackit(geo);

  tracker.trackdb()->display();

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Tracker.cc
//----------------------------------------------------------------------------//
