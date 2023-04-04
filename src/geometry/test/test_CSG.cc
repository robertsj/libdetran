//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_CSG.cc
 *  @brief Test of CSG classes
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "geometry/CSG.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>
#define COUT(c) std::cout << c << std::endl;

using namespace detran_geometry;
using namespace detran_utilities;
using std::cout;
using std::endl;

typedef QuadraticSurfaceFactory     QSF;
typedef Surface::vec_surface        vec_surface;
typedef CSG_Node::SP_node           SP_node;
typedef CSG_Node::vec_point         vec_point;
typedef std::vector<SP_node>        vec_node;

//----------------------------------------------------------------------------//
vec_surface get_surfaces()
{
  vec_surface surfaces(8);
  surfaces[0] = QSF::CreatePlaneX(0.0);
  surfaces[1] = QSF::CreatePlaneX(1.0);
  surfaces[2] = QSF::CreatePlaneY(0.0);
  surfaces[3] = QSF::CreatePlaneY(1.0);
  surfaces[4] = QSF::CreatePlaneZ(0.0);
  surfaces[5] = QSF::CreatePlaneZ(1.0);
  surfaces[6] = QSF::CreateSphere(0.5, 0.5, 0.5, 0.45);
  surfaces[7] = QSF::CreateSphere(1.3, 0.5, 0.5, 0.45);
  return surfaces;
}

//----------------------------------------------------------------------------//
vec_node get_primitives()
{
  // Surfaces and senses
  vec_surface surfaces = get_surfaces();
  vec_size_t sense(surfaces.size(), true);
  sense[1] = false;
  sense[3] = false;
  sense[5] = false;
  sense[6] = true;
  // Primitives
  vec_node primitives;
  for (int i = 0; i < 7; ++i)
  {
    SP_node tmp(new CSG_Primitive(surfaces[i], sense[i]));
    primitives.push_back(tmp);
  }
  return primitives;
}

//----------------------------------------------------------------------------//
SP_node get_box()
{
  vec_node primitives = get_primitives();
  SP_node box(new CSG_Intersection(primitives[0], primitives[1]));
  for (int i = 2; i < 6; ++i)
    box  = std::make_shared<CSG_Intersection>(box, primitives[i]);
  return box;
}

//----------------------------------------------------------------------------//
SP_node get_box_no_sphere()
{
  SP_node box = get_box();
  vec_node primitives = get_primitives();
  SP_node box_no_sphere(new CSG_Intersection(box, primitives[6]));
  return box_no_sphere;
}

//----------------------------------------------------------------------------//
TEST(CSG, Intersection)
{
  SP_node box = get_box();
  SP_node box_no_sphere = get_box_no_sphere();

  // Test contains
  EXPECT_TRUE(box->contains(Point(0.5, 0.5, 0.5)));;
  EXPECT_FALSE(box->contains(Point(2.0, 2.0, 2.0)));;
  EXPECT_FALSE(box_no_sphere->contains(Point(0.5, 0.5, 0.5)));;
  EXPECT_TRUE(box_no_sphere->contains(Point(0.01, 0.01, 0.01)));;

  // Test intersections
  Ray R(Point(-2.0, 0.5, 0.5), Point(1.0, 0.0, 0.0));

  vec_point points = box->intersections(R, 20.0);
  EXPECT_EQ(points.size(), 2);
  for (int i = 0; i < 2; ++i)
    COUT(points[i])
  EXPECT_NEAR(points[0].x(), 0.0, 1.0e-12);
  EXPECT_NEAR(points[1].x(), 1.0, 1.0e-12);

  points = box_no_sphere->intersections(R, 20.0);
  EXPECT_EQ(points.size(), 4);
  COUT(" ")
  for (int i = 0; i < 4; ++i)
    COUT(points[i])
  EXPECT_NEAR(points[0].x(), 0.00, 1.0e-12);
  EXPECT_NEAR(points[1].x(), 0.05, 1.0e-12);
  EXPECT_NEAR(points[2].x(), 0.95, 1.0e-12);
  EXPECT_NEAR(points[3].x(), 1.00, 1.0e-12);

}

//----------------------------------------------------------------------------//
TEST(CSG, Union)
{
  vec_surface surfaces = get_surfaces();
  SP_node sphere1(new CSG_Primitive(surfaces[6], false));
  SP_node sphere2(new CSG_Primitive(surfaces[7], false));
  SP_node double_sphere(new CSG_Union(sphere1, sphere2));

  // Test contains
  EXPECT_TRUE(double_sphere->contains(Point(0.7, 0.5, 0.5)));;
  EXPECT_TRUE(double_sphere->contains(Point(1.0, 0.5, 0.5)));;
  EXPECT_FALSE(double_sphere->contains(Point(1.9, 0.5, 0.5)));;

  // Test intersections
  Ray R(Point(-2.0, 0.5, 0.5), Point(1.0, 0.0, 0.0));
  vec_point points = double_sphere->intersections(R, 20.0);
  EXPECT_EQ(points.size(), 4);
  for (int i = 0; i < 4; ++i)
    COUT(points[i])
  EXPECT_NEAR(points[0].x(), 0.05, 1.0e-12);
  EXPECT_NEAR(points[1].x(), 0.85, 1.0e-12);
  EXPECT_NEAR(points[2].x(), 0.95, 1.0e-12);
  EXPECT_NEAR(points[3].x(), 1.75, 1.0e-12);
}

//----------------------------------------------------------------------------//
TEST(CSG, Difference)
{
  vec_node primitives = get_primitives();
  SP_node box = get_box();
  SP_node box_minus_sphere(new CSG_Difference(box, primitives[6]));

  // Test contains
  EXPECT_TRUE(box_minus_sphere->contains(Point(0.5, 0.5, 0.5)));;
  EXPECT_FALSE(box_minus_sphere->contains(Point(1.0, 0.5, 0.0)));;

  // Test intersections
  Ray R(Point(-2.0, 0.5, 0.5), Point(1.0, 0.0, 0.0));
  vec_point points = box_minus_sphere->intersections(R, 20.0);
  EXPECT_EQ(points.size(), 4);
  for (int i = 0; i < 4; ++i)
    COUT(points[i])
  EXPECT_NEAR(points[0].x(), 0.00, 1.0e-12);
  EXPECT_NEAR(points[1].x(), 0.05, 1.0e-12);
  EXPECT_NEAR(points[2].x(), 0.95, 1.0e-12);
  EXPECT_NEAR(points[3].x(), 1.00, 1.0e-12);
  EXPECT_FALSE(box_minus_sphere->contains(0.5 * (points[0] + points[1])));;
  EXPECT_TRUE(box_minus_sphere->contains(0.5 * (points[1] + points[2])));;
  EXPECT_FALSE(box_minus_sphere->contains(0.5 * (points[2] + points[3])));;
}

//----------------------------------------------------------------------------//
TEST(CSG, Translation)
{
  SP_node box = get_box();
  SP_node box_no_sphere = get_box_no_sphere();
  SP_node shifted(new CSG_Translation(box_no_sphere, Point(10, 0, 0)));

  // Test contains
  EXPECT_FALSE(box_no_sphere->contains(Point(0.50, 0.50, 0.50)));;
  EXPECT_TRUE(box_no_sphere->contains(Point(0.01, 0.01, 0.01)));;
  EXPECT_FALSE(shifted->contains(Point(10.50, 0.50, 0.50)));;
  EXPECT_TRUE(shifted->contains(Point(10.01, 0.01, 0.01)));;

  // Test intersections
  Ray R(Point(-2.0, 0.5, 0.5), Point(1.0, 0.0, 0.0));
  vec_point points = shifted->intersections(R, 20.0);
  EXPECT_EQ(points.size(), 4);
  for (int i = 0; i < points.size(); ++i)
    COUT(points[i])
  EXPECT_NEAR(points[0].x(), 10.00, 1.0e-12);
  EXPECT_NEAR(points[1].x(), 10.05, 1.0e-12);
  EXPECT_NEAR(points[2].x(), 10.95, 1.0e-12);
  EXPECT_NEAR(points[3].x(), 11.00, 1.0e-12);
}

//----------------------------------------------------------------------------//
//              end of test_Region.cc
//----------------------------------------------------------------------------//
