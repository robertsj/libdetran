//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_State.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  Test of State.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_State_basic)

// Detran headers
#include "TestDriver.hh"
#include "State.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"
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

int test_State_basic()
{

  // Test fixtures
  SP_mesh mesh          = mesh_2d_fixture();
  SP_quadrature quad    = quadruplerange_fixture();

  // Input
  State::SP_input input;
  input = new InputDB();
  input->put<string>("equation", "dd");
  input->put<int>("number_groups", 2);

  // State
  State::SP_state state;
  state = new State(input, mesh, quad);

  // Fill the fluxes.
  for (int cell = 0; cell < 1; cell++)
  {
    // Group 0
    (state->phi(0))[cell] = 1.23;
    // Group 1
    (state->phi(1))[cell] = 2.34;
  }

  State::moments_type phi_0 = state->phi(0);
  State::moments_type phi_1(state->phi(1));

  TEST(phi_0[0] == 1.23);
  TEST(phi_1[0] == 2.34);

  // Should not affect.
  phi_0[0] = 3.45;
  TEST((state->phi(0))[0] == 1.23);
  state->phi(0) = phi_0;
  TEST((state->phi(0))[0] == 3.45);

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_State.cc
//---------------------------------------------------------------------------//
