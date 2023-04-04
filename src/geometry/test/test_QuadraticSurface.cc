//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_QuadraticSurface.cc
 *  @brief Test of QuadraticSurface class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "geometry/QuadraticSurface.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/Constants.hh"
#include <cmath>

using namespace detran_geometry;
using namespace detran_utilities;
using std::cout;
using std::endl;

// Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J
TEST(QuadraticSurface, Plane)
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
    EXPECT_FALSE(x_plane.sense(r0));
    EXPECT_FALSE(y_plane.sense(r0));
    EXPECT_FALSE(y_plane.sense(r0));
    EXPECT_TRUE(x_plane.sense(r1));
    EXPECT_TRUE(y_plane.sense(r1));
    EXPECT_TRUE(y_plane.sense(r1));

    vec_point x_int = x_plane.intersections(Ray(r0, d));
    EXPECT_EQ(x_int.size(), 1);
    cout << x_int[0] << " " << x_plane.f(x_int[0]) << endl;
    EXPECT_NEAR(x_int[0].x(), 1.0, 1.0e-12);
    EXPECT_NEAR(x_int[0].y(), 1.0, 1.0e-12);
    EXPECT_NEAR(x_int[0].z(), 1.0, 1.0e-12);

    x_int = x_plane.intersections(Ray(r1, d));
    EXPECT_EQ(x_int.size(), 0);
  }
}

//----------------------------------------------------------------------------//
TEST(QuadraticSurface, Circle)
{
  typedef QuadraticSurfaceFactory QSF;
  QSF::SP_surface c = QSF::CreateCylinderZ(0.0, 0.0, 1.0);
  Point r0(0, 0.99, 0);
  Point r1(0, 1.01, 0);
  EXPECT_FALSE(c->sense(r0));
  EXPECT_TRUE(c->sense(r1));
  Point r(-2, -2, 0);
  Point d(1, 1, 0); d = d * (1.0 / distance(d));
  QuadraticSurface::vec_point points = c->intersections(Ray(r, d));
}

//----------------------------------------------------------------------------//
//              end of test_QuadraticSurface.cc
//----------------------------------------------------------------------------//
