//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_SiloOutput.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of PPMPlotter class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_PPMPlotter)

// Detran headers
#include "utilities/TestDriver.hh"
#include "ioutils/PPMPlotter.hh"
#include <iostream>

using namespace detran_test;
using namespace detran_utilities;
using namespace detran_ioutils;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_PPMPlotter(int argc, char *argv[])
{
  PPMPlotter plotter;
  int N = 1000;
  plotter.initialize(N, N, "test.ppm");
  for (int j = 0; j < N; ++j)
  {
    for (int i = 0; i < N; ++i)
    {
      double x = double(i)/double(N);
      double y = double(j)/double(N);
      double v = x;
      if (i < N/2 && j < N/2) v = -1.0;
      plotter.set_pixel(i, j, v);
    }
  }
  TEST(plotter.write());

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_SiloOutput.cc
//---------------------------------------------------------------------------//
