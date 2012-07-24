//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Point.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Point class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Point)

// Detran headers
#include "TestDriver.hh"
#include "Point.hh"

// Setup
/* ... */

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Point()
{
  Point p1(1.0, 2.0);
  TEST(soft_equiv(p1.x(), 1.0));
  TEST(soft_equiv(p1.y(), 2.0));
  Point p2 =  2.0 * p1;
  TEST(soft_equiv(p2.x(), 2.0));
  TEST(soft_equiv(p2.y(), 4.0));
  Point p3 = p1 + p2;
  TEST(soft_equiv(p3.x(), 3.0));
  TEST(soft_equiv(p3.y(), 6.0));
  Point p4 = p1 - p2;
  TEST(soft_equiv(p4.x(), -1.0));
  TEST(soft_equiv(p4.y(), -2.0));
  TEST(soft_equiv(detran::distance(p1, p2), 2.23606797749979));
  //
  p1 = Point(0.75, 0.00);
  p2 = Point(0.50, 0.50);
  TEST(soft_equiv(detran::distance(p2, p1), 0.559016994374947));
  TEST(soft_equiv(detran::distance(p1, p2), 0.559016994374947));


  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Point.cc
//---------------------------------------------------------------------------//
