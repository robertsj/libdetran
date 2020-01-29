//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_RegionFactory.cc
 *  @brief Test of RegionFactory struct
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                           \
        FUNC(test_RegionFactory_PinCell)    \
        FUNC(test_RegionFactory_Translated)

#include "TestDriver.hh"
#include "geometry/RegionFactory.hh"
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

int test_RegionFactory_PinCell(int argc, char *argv[])
{
  typedef RegionFactory RF;

  Point pitch(1.26);
  PinCell::vec_int mat_map(3, 0); mat_map[1] = 1; mat_map[2] = 2;
  PinCell::vec_dbl radii(2, 0.49);
  radii[1] = 0.54;
  PinCell::SP_pincell pin;
  Point center(0, 0);
  pin = PinCell::Create(pitch, mat_map, radii, PinCell::DIVISION_NONE, center);

  RF::vec_region regions;
  regions = RF::CreatePinCell(pin);


  return 0;
}

int test_RegionFactory_Translated(int argc, char *argv[])
{
  typedef RegionFactory RF;

  Point pitch(1.26);
  PinCell::vec_int mat_map(1, 0);
  PinCell::SP_pincell pin = PinCell::Create(pitch, mat_map);
  RF::vec_region regions = RF::CreatePinCell(pin);

  TEST( regions[0]->contains(Point(0.5, 0.5)));
  TEST(!regions[0]->contains(Point(1.5, 0.5)));

  RF::SP_region translated =
    RF::CreateTranslatedRegion(regions[0], Point(10, 10));

  TEST( translated->contains(Point(10.5, 10.5)));
  TEST(!translated->contains(Point(11.5, 10.5)));

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_RegionFactory.cc
//----------------------------------------------------------------------------//
