//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_LossOperator.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of LossOperator class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_LossOperator)

// Detran headers
#include "TestDriver.hh"
#include "LossOperator.hh"
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

// Test of basic public interface
int test_LossOperator()
{

  // Initialize PETSc
  int argc = 0;
  char **args;
  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);

  // Create input.
  LossOperator::SP_input inp(new InputDB());
  inp->put<string>("bc_left",  "reflect");
  inp->put<string>("bc_right", "reflect");

  // Create material.
  LossOperator::SP_material mat(new Material(2, 1, false));
  mat->set_sigma_t(0, 0, 0.1890);
  mat->set_sigma_t(0, 1, 1.4633);
  mat->set_sigma_s(0, 0, 0, 0.1507);
  mat->set_sigma_s(0, 0, 1, 0.0000);
  mat->set_sigma_s(0, 1, 0, 0.0380);
  mat->set_sigma_s(0, 1, 1, 1.4536);
  mat->set_diff_coef(0, 0, 0.1890/3.0);
  mat->set_diff_coef(0, 1, 1.4633/3.0);
  mat->finalize();

  // phi0 = 1 / (0.1890-0.1507)
  // phi1 = (1 + phi0 * 0.0380) / (1.4633-1.4536)

  // Create mesh.
  vec_dbl cm(2, 0.0);
  cm[1] = 10.0;
  vec_int fm(1, 10);
  vec_int mt(1, 0);
  Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));
  //Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  //Mesh3D::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create the operator.
  LossOperator M(inp, mat, mesh);

  // Create a right hand side.
  Vec x;
  Vec y;
  int n = mat->number_groups() * mesh->number_cells();
  VecCreateSeq(PETSC_COMM_SELF, n, &x);
  VecDuplicate(x, &y);
  VecSet(x, 1.0);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  M.multiply(x, y);
  M.display();
//
//  double *y_a;
//  VecGetArray(y, &y_a);
//  TEST(soft_equiv(y_a[0], 0.385714285714286));
//  for (int i = 1; i < 9; i++)
//    TEST(soft_equiv(y_a[i], 0.1000000000000000));
//  TEST(soft_equiv(y_a[9], 0.385714285714286));
//  VecRestoreArray(y, &y_a);
 // PetscFinalize();

  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp, M.get_operator(), M.get_operator(), SAME_NONZERO_PATTERN);
  KSPSolve(ksp, x, y);
  VecView(y, PETSC_VIEWER_STDOUT_SELF);
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_LossOperator.cc
//---------------------------------------------------------------------------//
