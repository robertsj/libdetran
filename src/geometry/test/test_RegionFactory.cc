//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_RegionFactory.cc
 *  @brief Test of RegionFactory struct
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "geometry/RegionFactory.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_geometry;
using namespace detran_utilities;
using std::cout;
using std::endl;

TEST(RegionFactory, PinCell)
{
  typedef RegionFactory RF;

  Point pitch(1.26);
  PinCell::vec_int mat_map(3, 0); mat_map[1] = 1; mat_map[2] = 2;
  PinCell::vec_dbl radii(2, 0.49);
  radii[1] = 0.54;
  PinCell::SP_pincell pin;
  Point center(0, 0);
  pin = std::make_shared<PinCell>(pitch, mat_map, radii, PinCell::DIVISION_NONE, center);

  RF::vec_region regions;
  regions = RF::CreatePinCell(pin);

}

TEST(RegionFactory, Translated)
{
  typedef RegionFactory RF;

  Point pitch(1.26);
  PinCell::vec_int mat_map(1, 0);
  PinCell::SP_pincell pin = std::make_shared<PinCell>(pitch, mat_map);
  RF::vec_region regions = RF::CreatePinCell(pin);

  EXPECT_TRUE(regions[0]->contains(Point(0.5, 0.5)));;
  EXPECT_FALSE(regions[0]->contains(Point(1.5, 0.5)));;

  RF::SP_region translated =
    RF::CreateTranslatedRegion(regions[0], Point(10, 10));

  EXPECT_TRUE(translated->contains(Point(10.5, 10.5)));;
  EXPECT_FALSE(translated->contains(Point(11.5, 10.5)));;
}

//----------------------------------------------------------------------------//
//              end of test_RegionFactory.cc
//----------------------------------------------------------------------------//
