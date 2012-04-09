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
#include "Sweeper.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "angle/test/quadrature_fixture.hh"

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

int test_Sweeper2D_basic()
{
  typedef Sweeper<_2D> Sweeper2D;

  // Test fixtures
  SP_mesh mesh          = mesh_2d_fixture();
  SP_material mat       = material_fixture_1g();
  SP_quadrature quad    = quadruplerange_fixture();

  // Input
  Sweeper2D::SP_input input;
  input = new InputDB();
  input->put<string>("equation", "dd");
  input->put<int>("number_groups", 2);

  // State
  Sweeper2D::SP_state state;
  state = new State(input, mesh, quad);

  // Boundary
  Boundary<_2D>::SP_boundary bound;
  bound = new Boundary<_2D>(input, mesh, quad);

  // Sweeper
//  Sweeper2D sweeper(input, mesh, mat,
//                    quad, state, bound, source);
//
//  // Get moments and make source.
//  State::moments_type phi = state->phi(0);
//  State::moments_type source(phi);
//  source.assign(source.size(), 1.0);
//
//  // Sweep.
//  sweeper.setup_group(0);
  //sweeper.sweep(phi);

  // Tests.

  return 0;
}
