//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Sweeper2D.cc
 *  @brief test_Sweeper2D
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Sweeper2D_basic)

#include "utilities/TestDriver.hh"
#include "Sweeper2D.hh"
#include "Equation_DD_2D.hh"
#include "Equation_SD_2D.hh"
#include "Equation_SC_2D.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "angle/QuadratureFactory.hh"

using namespace detran;
using namespace detran_angle;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_Sweeper2D_basic(int argc, char *argv[])
{
  typedef Sweeper2D<Equation_DD_2D> Sweeper_T;

  // Test fixtures
  SP_mesh mesh          = mesh_2d_fixture();
  SP_material mat       = material_fixture_1g();

  // Input
  Sweeper_T::SP_input input;
  input = new InputDB();
  input->put<std::string>("equation", "dd");
  input->put<int>("number_groups", 2);

  QuadratureFactory::SP_quadrature quad = QuadratureFactory::build(input, 2);

  // State
  Sweeper_T::SP_state state;
  state = new State(input, mesh, quad);

  // Boundary
  Sweeper_T::SP_boundary bound;
  bound = new Sweeper_T::Boundary_T(input, mesh, quad);

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

//----------------------------------------------------------------------------//
//              end of test_Sweeper2D.cc
//----------------------------------------------------------------------------//
