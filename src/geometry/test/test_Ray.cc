//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_Ray.cc
 *  @author Jeremy Roberts
 *  @note   Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_Ray)

#include "TestDriver.hh"
#include "Ray.hh"

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

int test_Ray(int argc, char *argv[])
{
  Point origin(-2.0, 0.5, 0.5);
  Point direction(1.0, 0.0, 0.0);
  Ray r(origin, direction);
  double tt = 2.0;
  Point p =  r.origin + tt*r.direction;
  cout << origin << endl;
  cout << p      << endl;

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Ray.cc
//---------------------------------------------------------------------------//
