//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_BoundaryMOC.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of BoundaryMOC
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_BoundaryMOC_2D)

// Detran headers
#include "TestDriver.hh"
#include "BoundaryMOC.hh"
#include "Mesh2D.hh"
#include "Uniform.hh"
#include "Tracker.hh"

// Setup
//

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

int test_BoundaryMOC_2D()
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh0(new Mesh2D(fm, fm, cm, cm, mt));

  // Create quadrature.  1 azimuth, 3 space, 1 polar.
  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));
  quad->display_tracks();
  // Create tracker
  Tracker tracker(mesh0, quad);
  Tracker::SP_trackdb tracks = tracker.trackdb();

  // Get the tracked mesh.
  Mesh::SP_mesh mesh = tracker.meshmoc();

  TEST(mesh);
  TEST(quad);


  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Equation_SC_MOC.cc
//---------------------------------------------------------------------------//
