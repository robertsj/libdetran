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

int test_CSG(int argc, char *argv[])
{
  typedef QuadraticSurfaceFactory QSF;
  typedef CSG_Node::SP_node  SP_node;

  // Create a box via surfaces (simulates Region innards)
  Surface::vec_surface surfaces(7);

  // Surfaces
  surfaces[0] = QSF::CreatePlaneX(0.0);
  surfaces[1] = QSF::CreatePlaneX(1.0);
  surfaces[2] = QSF::CreatePlaneY(0.0);
  surfaces[3] = QSF::CreatePlaneY(1.0);
  surfaces[4] = QSF::CreatePlaneZ(0.0);
  surfaces[5] = QSF::CreatePlaneZ(1.0);
  surfaces[6] = QSF::CreateSphere(0.5, 0.5, 0.5, 0.45);

  // Senses
  vec_size_t sense(surfaces.size(), 1);
  sense[1] = false;
  sense[3] = false;
  sense[5] = false;
  sense[6] = true;

  // Primitives
  std::vector<SP_node> primitives;
  for (int i = 0; i < 7; ++i)
  {
    SP_node tmp(new CSG_Primitive(surfaces[i], sense[i]));
    primitives.push_back(tmp);
  }

  // Box node
  SP_node box(new CSG_Intersection(primitives[0], primitives[1]));
  for (int i = 2; i < 6; ++i)
  {
    box = new CSG_Intersection(box, primitives[i]);
  }
  TEST(box->contains(Point(0.5, 0.5, 0.5)));
  TEST(!box->contains(Point(2.0, 2.0, 2.0)));

  // In box, outside sphere
  SP_node box_no_sphere(new CSG_Intersection(box, primitives[6]));
  TEST(!box_no_sphere->contains(Point(0.5, 0.5, 0.5)));
  TEST(box_no_sphere->contains(Point(0.01, 0.01, 0.01)));

  // Shift the whole thing (0, 0, 0) --> (10, 10, 10)
  SP_node shifted(new CSG_Translation(box_no_sphere, Point(10,10,10)));
  TEST(!shifted->contains(Point(10.5, 10.5, 10.5)));
  TEST(shifted->contains(Point(10.01, 10.01, 10.01)));

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Region.cc
//----------------------------------------------------------------------------//
