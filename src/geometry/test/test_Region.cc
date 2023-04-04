//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Region.cc
 *  @brief Test of Region class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "geometry/Region.hh"
#include "geometry/CSG.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_geometry;
using namespace detran_utilities;
using std::cout;
using std::endl;

TEST(Region, Basic)
{
  typedef QuadraticSurfaceFactory QSF;

  Surface::SP_surface W = QSF::CreatePlaneX(0.0);
  Surface::SP_surface E = QSF::CreatePlaneX(1.0);
  Surface::SP_surface S = QSF::CreatePlaneY(0.0);
  Surface::SP_surface N = QSF::CreatePlaneY(1.0);
  Surface::SP_surface B = QSF::CreatePlaneZ(0.0);
  Surface::SP_surface T = QSF::CreatePlaneZ(1.0);
  Surface::SP_surface C = QSF::CreateCylinderZ(0.5, 0.5, 0.4);

  Region::SP_region region0 =
   std::make_shared<Region>(0, Point(0, 0, 0) - 0.1, Point(1, 1, 1) + 0.1);
  region0->append(C, false);

  Region::SP_region region1 =
   std::make_shared<Region>(1, Point(0, 0, 0) - 0.1, Point(1, 1, 1) + 0.1);
  region1->append(W, true);
  region1->append(E, false);
  region1->append(S, true);
  region1->append(N, false);
  region1->append(B, true);
  region1->append(T, false);
  region1->append(C, true);

  // test contains
  EXPECT_TRUE(region0->contains(Point(0.5, 0.5, 0.5)));;
  EXPECT_TRUE(region1->contains(Point(0.1, 0.1, 0.1)));;

  // test bounding box intersections
  Ray r(Point(-2.0, 0.5, 0.5), Point(1, 0, 0));
  EXPECT_TRUE(region0->intersects_bounding_box(r, 20.0));;
}

//----------------------------------------------------------------------------//
//              end of test_Region.cc
//----------------------------------------------------------------------------//
