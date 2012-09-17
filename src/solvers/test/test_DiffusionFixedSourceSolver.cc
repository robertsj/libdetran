//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_DiffusionFixedSourceSolver.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  test_DiffusionFixedSourceSolver class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                                \
        FUNC(test_DiffusionFixedSourceSolver_1D)

// Detran headers
#include "TestDriver.hh"
#include "DiffusionFixedSourceSolver.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"
#include "external_source/ConstantSource.hh"

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "external_source/test/external_source_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  RUN(argc, argv);
  PetscFinalize();
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_DiffusionFixedSourceSolver_1D(int argc, char *argv[])
{

  // Material (kinf = 5/4 = 1.25)
  SP_material mat(new Material(1, 1));
  mat->set_sigma_t(0, 0,    1.0);
  mat->set_sigma_s(0, 0, 0, 0.6);
  mat->set_sigma_f(0, 0,    0.5);
  mat->set_diff_coef(0, 0,  0.33);
  mat->finalize();
  // Mesh
  vec_dbl cm(2, 0.0); cm[1] = 10.0;
  vec_int fm(1, 10);
  vec_int mat_map(1, 0);
  SP_mesh mesh(new Mesh1D(fm, cm, mat_map));

  // Constant unit source.
  ConstantSource::SP_externalsource q_e(new ConstantSource(1, mesh, 1.0));

  // Input
  DiffusionFixedSourceSolver<_1D>::SP_input input;
  input = new InputDB();
  input->put<int>(     "number_groups",         1);
  input->put<string>(  "bc_west",               "reflect");
  input->put<string>(  "diffusion_fixed_type",  "reflect");

  // State
  DiffusionFixedSourceSolver<_1D>::SP_state state(new State(input, mesh));

  // solver
  cout << " Building solver... " << endl;
  DiffusionFixedSourceSolver<_1D> solver(input, mat, mesh, state);

  // Build source
  cout << " Building source... " << endl;
  solver.build_source(q_e);

  // Solve.
  cout << " Solving... " << endl;
  solver.solve();

  for (int i = 0; i < 10; i++)
    cout << " i = " << i << " phi = " << state->phi(0)[i] << endl;
  DiffusionFixedSourceSolver<_1D>::SP_boundary boundary;
  boundary = solver.boundary();

//  BoundaryValue<D>::value
//              (J(surface, g, Boundary_T::OUT), ijk[dim1], ijk[dim2]) =
//                ((2.0 * DC) * (*d_phi)[row] + (W - 4.0*DC) * Jinc) / (4.0 * DC + W);

  cout << (*boundary)(0, 0, 1) << endl;
  cout << (*boundary)(0, 0, 0) << endl;
  cout << (*boundary)(1, 0, 1) << endl;
  cout << (*boundary)(1, 0, 0) << endl;

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_DiffusionFixedSourceSolver.cc
//---------------------------------------------------------------------------//
