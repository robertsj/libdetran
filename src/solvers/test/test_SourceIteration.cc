//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_SourceIteration.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  test_SourceIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                      \
        FUNC(test_SourceIteration_2D   )

// Detran headers
#include "TestDriver.hh"
#include "SourceIteration.hh"

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "transport/test/external_source_fixture.hh"

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

int test_SourceIteration_2D()
{

  // Test fixtures
  SP_material mat       = material_fixture_2g();
  SP_mesh mesh          = mesh_2d_fixture();
  SP_quadrature quad    = quadruplerange_fixture();
  SP_source q_e         = constant_source_fixture(2);
  SourceIteration<_2D>::SP_fissionsource q_f;

  // Input
  SourceIteration<_2D>::SP_input input;
  input = new InputDB();
  input->put<string>("equation", "dd");
  input->put<int>("number_groups", 2);

  // State
  SourceIteration<_2D>::SP_state state;
  state = new State(input, mesh, quad);

  // Boundary
  SourceIteration<_2D>::SP_boundary bound;
  bound = new Boundary<_2D>(input, mesh, quad);

  // SI
  SourceIteration<_2D> solver(input, state, mesh, mat,
                              quad, bound, q_e, q_f);


  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_SourceIteration.cc
//---------------------------------------------------------------------------//
