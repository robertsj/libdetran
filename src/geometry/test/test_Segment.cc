//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Segment.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Segment class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Segment)

// Detran headers
#include "TestDriver.hh"
#include "Segment.hh"

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

int test_Segment()
{
  Segment s(0, 1.0);
  TEST(s.region() == 0);
  TEST(soft_equiv(s.length(), 1.0));
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Segment.cc
//---------------------------------------------------------------------------//
