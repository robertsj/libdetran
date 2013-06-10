//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Region.cc
 *  @brief Test of Region class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                          \
        FUNC(test_Region)

#include "TestDriver.hh"
#include "geometry/Region.hh"
#include "geometry/CSG.hh"
#include <cmath>

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Given a ray and solid, ray cast returns the segment lengths and start
// points
void ray_cast(Region::SP_region region,
              const Point      &start,
              const Point      &end)
{

}

//
CSG_Node::vec_point ray_cast(Region::SP_node  node,
                             const Point      &start,
                             const Point      &end)
{
  // define direction
  double length = distance(end, start);
  Point d = (end - start) / length;
  // get intersections
  CSG_Node::vec_point points = node->intersections(start, d, length);
  return points;
}


int test_Region(int argc, char *argv[])
{
  // Create a box via surfaces

  // Surfaces
  Surface::SP_surface W(new PlaneX(0.0));
  Surface::SP_surface E(new PlaneX(1.0));
  Surface::SP_surface S(new PlaneY(0.0));
  Surface::SP_surface N(new PlaneY(0.0));
  Surface::SP_surface B(new PlaneZ(0.0));
  Surface::SP_surface T(new PlaneZ(1.0));

  // Node
  CSG_Node::SP_node n_W(new CSG_Primitive(W, true));
  CSG_Node::SP_node n_E(new CSG_Primitive(E, false));
  CSG_Node::SP_node n_S(new CSG_Primitive(S, true));
  CSG_Node::SP_node n_N(new CSG_Primitive(N, false));
  CSG_Node::SP_node n_B(new CSG_Primitive(B, true));
  CSG_Node::SP_node n_T(new CSG_Primitive(T, false));

  CSG_Node::SP_node node(new Union())

  // Region [W ^ -E ^ S ^ -N ^ B ^ -T]
  // Region box((W, Intersect

  return 0;

}

//----------------------------------------------------------------------------//
//              end of test_Region.cc
//----------------------------------------------------------------------------//
