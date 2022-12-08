//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_InnerGMRES.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  test_InnerGMRES class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                      \
        FUNC(test_InnerGMRES_1D)

// Detran headers
#include "TestDriver.hh"
#include "InnerGMRES.hh"
#include "Mesh1D.hh"
#include "external_source/ConstantSource.hh"

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "external_source/test/external_source_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_angle;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_InnerGMRES_1D_actual();

int test_InnerGMRES_1D(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Call actual test.
  int result = test_InnerGMRES_1D_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;
}

int test_InnerGMRES_1D_actual()
{

  // Test fixtures.  Material 0 has sigma = 1, c = 0.9.
  SP_material mat = material_fixture_1g();
  mat->set_sigma_s(0, 0, 0, 0.99999);
  // Mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 10.0;
  vec_int fm(1, 100);
  vec_int mat_map(1, 0);
  SP_mesh mesh;
  mesh  = std::make_shared<Mesh1D>(fm, cm, mat_map);

  // Quadrature
  SP_quadrature quad = gausslegendre_fixture();

  // Constant unit source.
  ConstantSource::SP_externalsource
    q_e(new ConstantSource(1, mesh, 1.0, quad));

  // Input
  InnerGMRES<_1D>::SP_input input;
  input  = std::make_shared<InputDB>();
  input->put<string>(  "equation",          "dd");
  input->put<int>(     "number_groups",     1);
  input->put<int>(     "inner_max_iters",   100);
  input->put<double>(  "inner_tolerance",   1e-8);
  input->put<int>(     "inner_use_pc",      1);
  input->put<string>(     "bc_left",      "reflect");

  // State
  InnerGMRES<_1D>::SP_state state(new State(input, mesh, quad));

  // Fission source (uninitialized)
  InnerGMRES<_1D>::SP_fissionsource q_f(new FissionSource(state, mesh, mat));

  // Boundary
  InnerGMRES<_1D>::SP_boundary bound(new BoundarySN<_1D>(input, mesh, quad));

  // SI
  InnerGMRES<_1D> solver(input, state, mesh, mat, quad, bound, q_e, q_f);

  // Solve.
  solver.solve(0);

  for (int i = 0; i < 10; i++)
    cout << " i = " << i << " phi = " << state->phi(0)[i] << endl;

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_InnerGMRES.cc
//---------------------------------------------------------------------------//
