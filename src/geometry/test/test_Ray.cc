//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_Ray.cc
 *  @author Jeremy Roberts
 *  @note   Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Ray.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using std::cout;
using std::endl;

TEST(Ray, Basic)
{
  Point origin(-2.0, 0.5, 0.5);
  Point direction(1.0, 0.0, 0.0);
  Ray r(origin, direction);
  double tt = 2.0;
  Point p =  r.origin + tt*r.direction;
  cout << origin << endl;
  cout << p      << endl;
}

//---------------------------------------------------------------------------//
//              end of test_Ray.cc
//---------------------------------------------------------------------------//
