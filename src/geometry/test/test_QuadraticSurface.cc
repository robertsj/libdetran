//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_QuadraticSurface.cc
 *  @brief Test of QuadraticSurface class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                          \
        FUNC(test_QuadraticSurface_plane)  \
        FUNC(test_QuadraticSurface_circle)

#include "TestDriver.hh"
#include "geometry/QuadraticSurface.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/Constants.hh"
#include <cmath>

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J
int test_QuadraticSurface_plane(int argc, char *argv[])
{

  typedef Surface::vec_point vec_point;

  //--------------------------------------------------------------------------//
  // Test of x, y, and z planes, each at 1.0
  //--------------------------------------------------------------------------//

  {

    QuadraticSurface x_plane(0, 0, 0,   0, 0, 0,   1, 0, 0,  -1.0);
    QuadraticSurface y_plane(0, 0, 0,   0, 0, 0,   0, 1, 0,  -1.0);
    QuadraticSurface z_plane(0, 0, 0,   0, 0, 0,   0, 0, 1,  -1.0);

    // The starting point is at the origin, and we're going in a 45
    // degree angle w/r to each axis.
    Point r0(0, 0, 0);
    Point r1(2, 3, 4);
    Point d(1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0));

    cout << x_plane.f(r0) << endl;
    cout << x_plane.f(r1) << endl;
    cout << distance(d)   << endl;
    cout << d  << endl;
    cout << r0 << endl;
    TEST(x_plane.sense(r0) == false);
    TEST(y_plane.sense(r0) == false);
    TEST(y_plane.sense(r0) == false);
    TEST(x_plane.sense(r1) == true);
    TEST(y_plane.sense(r1) == true);
    TEST(y_plane.sense(r1) == true);

    vec_point x_int = x_plane.intersections(Ray(r0, d));
    TEST(x_int.size() == 1);
    cout << x_int[0] << " " << x_plane.f(x_int[0]) << endl;
    TEST(soft_equiv(x_int[0].x(), 1.0));
    TEST(soft_equiv(x_int[0].y(), 1.0));
    TEST(soft_equiv(x_int[0].z(), 1.0));

    x_int = x_plane.intersections(Ray(r1, d));
    TEST(x_int.size() == 0);
  }

  return 0;
}

//----------------------------------------------------------------------------//
int test_QuadraticSurface_circle(int argc, char *argv[])
{
  typedef QuadraticSurfaceFactory QSF;
  QSF::SP_surface c = QSF::CreateCylinderZ(0.0, 0.0, 1.0);
  Point r0(0, 0.99, 0);
  Point r1(0, 1.01, 0);
  TEST(c->sense(r0) == false);
  TEST(c->sense(r1) == true);
  Point r(-2, -2, 0);
  Point d(1, 1, 0); d = d * (1.0 / distance(d));
  QuadraticSurface::vec_point points = c->intersections(Ray(r, d));
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_QuadraticSurface.cc
//----------------------------------------------------------------------------//
