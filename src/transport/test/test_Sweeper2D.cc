//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Sweeper2D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of test_Sweeper2D
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Sweeper2D_basic)

// Detran headers
#include "TestDriver.hh"
#include "Sweeper2D.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "angle/test/quadrature_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace detran_utils;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Sweeper2D_basic()
{

  // Test fixtures
  SP_mesh mesh      = mesh_2d_fixture();
  SP_material mat   = material_fixture_1g();
  SP_quadrature q   = quadruplerange_fixture();

  // Input
  Sweeper2D::SP_input input;
  input = new InputDB();
  input->put<string>("equation", "dd");
  input->put<int>("number_groups", 2);

  // State
  Sweeper2D::SP_state state;
  state = new State(input, mesh, q);

  // Sweeper
  Sweeper2D sweeper(input, mesh, mat, q, state);

  // Get moments and make source.
  State::moments_type phi = state->phi(0);
  State::moments_type source(phi);
  source.assign(source.size(), 1.0);

  // Sweep.
  sweeper.sweep(phi, source);

  // Tests.

  return 0;
}
