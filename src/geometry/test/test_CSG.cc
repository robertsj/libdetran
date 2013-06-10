//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_CSG.cc
 *  @brief Test of CSG classes
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                          \
        FUNC(test_CSG)

#include "TestDriver.hh"
#include "geometry/CSG.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "callow/utils/Initialization.hh"
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

// Get segments
CSG_Node::vec_point ray_cast(CSG_Node::SP_node  node,
                             const Point       &start,
                             const Point       &end)
{
  // define direction
  double length = distance(end, start);
  Point d = (end - start) * (1.0/length);

  // get all intersections with the node constituents
  CSG_Node::vec_point points = node->intersections(start, d, length);

  // using the midpoint of each consecutive pair, determine the midpoint
  // and use this to check if it's in or out of the node.  This is a
  // simple way to handle concave regions, multiple adjoining body regions, etc.

  CSG_Node::vec_point midpoints;
  vec_dbl segments;
  for (int i = 0; i < points.size() - 1; ++i)
  {
    // make sure we have proper order
    Assert(distance(points[i], start) < distance(points[i+1], start));
    Point mid = 0.5*(points[i]+points[i+1]);
    std::cout << "mid=" << mid << std::endl;
    if (node->contains(mid))
    {
      segments.push_back(distance(points[i], points[i+1]));
      midpoints.push_back(mid);
    }
  }
  return midpoints;
}


int test_CSG(int argc, char *argv[])
{
  typedef QuadraticSurfaceFactory QSF;
  // Create a box via surfaces (simulates Region innards)
  Surface::vec_surface surfaces(6);

  // Surfaces
  surfaces[0] = QSF::PlaneX(0.0);
  surfaces[1] = QSF::PlaneX(1.0);
  surfaces[2] = QSF::PlaneY(0.0);
  surfaces[3] = QSF::PlaneY(1.0);
  surfaces[4] = QSF::PlaneZ(0.0);
  surfaces[5] = QSF::PlaneZ(1.0);

  // Senses
  vec_size_t sense(surfaces.size(), 1);
  sense[1] = false;
  sense[3] = false;
  sense[5] = false;

  // Primitives
  std::vector<CSG_Node::SP_node> primitives;
  for (int i = 0; i < surfaces.size(); ++i)
  {
    CSG_Node::SP_node tmp(new CSG_Primitive(surfaces[i], sense[i]));
    primitives.push_back(tmp);
  }

  // Box node
  CSG_Node::SP_node node(new CSG_Intersection(primitives[0], primitives[1]));
  for (int i = 2; i < surfaces.size(); ++i)
  {
    node = new CSG_Intersection(node, primitives[i]);
  }

  TEST(node->contains(Point(0.5, 0.5, 0.5)));
  TEST(!node->contains(Point(2.0, 2.0, 2.0)));

  CSG_Node::vec_point points =
    ray_cast(node, Point(-0.5, 0.5, 0.0), Point(1.5, 0.5, 0.0));
  std::cout << "final midpoints = " << std::endl;
  for (int i = 0; i < points.size(); ++i)
  {
    std::cout << points[i] << std::endl;
  }
  return 0;

}

//----------------------------------------------------------------------------//
//              end of test_Region.cc
//----------------------------------------------------------------------------//
