//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_DiffusionEigensolver.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of DiffusionEigensolver class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_DiffusionEigensolver)

// Detran headers
#include "TestDriver.hh"
#include "DiffusionEigensolver.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"

#include <iostream>

// Setup
//#include "material_fixture.hh"
using namespace detran;
using namespace detran_diffusion;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_DiffusionEigensolver_actual();

// Test of basic public interface
int test_DiffusionEigensolver(int argc, char *argv[])
{

  // Initialize SLEPc
  SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  cout << " number args " << argc << endl;

  // Run the actual test.
  int result = test_DiffusionEigensolver_actual();

  // Finalize SLEPc
  SlepcFinalize();

  return result;

}

int test_DiffusionEigensolver_actual()
{

  // Create input.
  DiffusionEigensolver::SP_input inp(new InputDB());
  inp->put<string>("bc_left",           "vacuum");
  inp->put<string>("bc_right",          "reflect");
  inp->put<int>("number_groups",        2);
  inp->put<int>("diffusion_max_iters",  100);

  // Create material.
  DiffusionEigensolver::SP_material mat(new Material(2, 1, false));
  mat->set_sigma_t(0, 0, 0.2263);
  mat->set_sigma_t(0, 1, 1.0119);
  mat->set_sigma_f(0, 0, 0.0067);
  mat->set_sigma_f(0, 1, 0.1241);
  mat->set_chi(0, 0, 1.0);
  mat->set_chi(0, 1, 0.0);
  mat->set_sigma_s(0, 0, 0, 0.2006);
  mat->set_sigma_s(0, 0, 1, 0.0000);
  mat->set_sigma_s(0, 1, 0, 0.0161);
  mat->set_sigma_s(0, 1, 1, 0.9355);
  mat->set_diff_coef(0, 0, 0.2263 / 3.0);
  mat->set_diff_coef(0, 1, 1.0119 / 3.0);
  mat->finalize();

  // Create mesh.
  vec_dbl cm(2, 0.0);
  cm[1] = 10.0;
  vec_int fm(1, 10);
  vec_int mt(1, 0);
  Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));
  //Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  //Mesh3D::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create state.
  DiffusionEigensolver::SP_state state(new State(inp, mesh));

  // Create the operator.
  DiffusionEigensolver solver(inp, mat, mesh, state);

  // Solve
  solver.solve();

  for (int i = 0; i < mesh->number_cells(); i++)
  {
    cout << state->phi(0)[i] << endl;
  }

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_DiffusionEigensolver.cc
//---------------------------------------------------------------------------//
