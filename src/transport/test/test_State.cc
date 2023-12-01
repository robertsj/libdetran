//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_State.cc
 *  @brief Test of State.
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "State.hh"
#include "geometry/test/mesh_fixture.hh"
#include "angle/QuadratureFactory.hh"

using namespace detran;
using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;

TEST(State, Basic)
{

  SP_mesh mesh          = mesh_2d_fixture();

  // Input
  State::SP_input input;
  input  = std::make_shared<InputDB>();
  input->put<std::string>("equation", "dd");
  input->put<int>("number_groups", 2);

  QuadratureFactory::SP_quadrature quad = QuadratureFactory::build(input, 2);

  // State
  State::SP_state state;
  state  = std::make_shared<State>(input, mesh, quad);

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

  EXPECT_EQ(phi_0[0], 1.23);
  EXPECT_EQ(phi_1[0], 2.34);

  // Should not affect.
  phi_0[0] = 3.45;
  EXPECT_EQ((state->phi(0))[0], 1.23);
  state->phi(0) = phi_0;
  EXPECT_EQ((state->phi(0))[0], 3.45);
}

//----------------------------------------------------------------------------//
//              end of test_State.cc
//----------------------------------------------------------------------------//
