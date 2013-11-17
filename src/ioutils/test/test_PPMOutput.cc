//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_PPMOutput.cc
 *  @brief Test of PPMPOutput class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_PPMOutput_mesh)  \
        FUNC(test_PPMOutput_geo)   \

#include "utilities/TestDriver.hh"
#include "ioutils/PPMOutput.hh"
#include "ioutils/ColorMap.hh"
#include "geometry/test/csg_fixture.hh"
#include "callow/utils/Initialization.hh"
#include <iostream>

using namespace detran_test;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_ioutils;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_PPMOutput_mesh(int argc, char *argv[])
{

  return 0;
}

int test_PPMOutput_geo(int argc, char *argv[])
{
//  // simple pin cell
//  {
//    PPMOutput::SP_geometry geo = test_2D_pincell_simple();
//    PPMOutput ppm("test_2D_pincell_simple");
//    ppm.initialize(geo, 0.01);
//    ppm.draw_geometry(geo, true);
//  }
//
//  // complex pin cell
//  {
//    PPMOutput::SP_geometry geo = test_2D_pincell_complex();
//    PPMOutput ppm("test_2D_pincell_complex");
//    ppm.initialize(geo, 0.001);
//    ppm.draw_geometry(geo, true);
//  }

  // complex pin cell
  {
    PPMOutput::SP_geometry geo = test_2D_pincell_via_factory(3, 3);
    //PPMOutput::SP_geometry geo = test_2D_assembly(8);
    PPMOutput ppm("test_2D_pincell_via_factory");
    ppm.initialize(geo, 0.01);
    ppm.draw_geometry(geo, true, ColorMap::RANDOM);
  }


  return 0;
}

//---------------------------------------------------------------------------//
//              end of PPMPOutput.cc
//---------------------------------------------------------------------------//
