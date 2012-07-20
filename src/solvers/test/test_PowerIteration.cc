//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_PowerIteration.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  test_PowerIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                      \
        FUNC(test_PowerIteration_2D   )

// Detran headers
#include "TestDriver.hh"
#include "PowerIteration.hh"

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "transport/test/external_source_fixture.hh"

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

int test_PowerIteration_2D()
{

  // Test fixtures
  SP_material mat = material_fixture_1g();

  // Mesh
  detran::vec_dbl cm(3, 0.0);
  cm[1] = 1.0;
  cm[2] = 2.0;
  detran::vec_int fm(2, 5);
  detran::vec_int mat_map(4, 2);
  SP_mesh mesh;
  mesh = new detran::Mesh2D(fm, fm, cm, cm, mat_map);

  // Quadrature
  SP_quadrature quad = quadruplerange_fixture();

  // Empty source.
  ExternalSource::SP_source q_e;

  // Input
  PowerIteration<_2D>::SP_input input;
  input = new InputDB();
  input->put<string>(  "equation",          "sc");
  input->put<int>(     "number_groups",     1);
  input->put<int>(     "inner_max_iters",   100);
  input->put<double>(  "inner_tolerance",   1e-12);
  input->put<int>(     "inner_print_out",   0);
  input->put<int>(     "outer_max_iters",   0);
  input->put<double>(  "outer_tolerance",   1e-12);
  input->put<int>(     "outer_print_out",   0);
  input->put<int>(     "eigen_max_iters",   100);
  input->put<double>(  "eigen_tolerance",   1e-14);
  input->put<string>(  "bc_left",           "reflect");
  input->put<string>(  "bc_right",          "reflect");
  input->put<string>(  "bc_bottom",         "reflect");
  input->put<string>(  "bc_top",            "reflect");

  // State
  PowerIteration<_2D>::SP_state state;
  state = new State(input, mesh, quad);

  // Fission source
  PowerIteration<_2D>::SP_fissionsource q_f;
  q_f = new FissionSource(state, mesh, mat);

  // Boundary
  PowerIteration<_2D>::SP_boundary bound;
  bound = new Boundary<_2D>(input, mesh, quad);

  // SI
  PowerIteration<_2D> solver(input, state, mesh, mat,
                             quad, bound, q_f);

  // Solve.
  solver.solve();
  TEST(soft_equiv(state->eigenvalue(), 1.0));
  TEST(soft_equiv(state->phi(0)[0], state->phi(0)[1]));

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_PowerIteration.cc
//---------------------------------------------------------------------------//
