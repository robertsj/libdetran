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

#include "utilities/TestDriver.hh"
#include "diffusion/LossOperator.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include <iostream>

// Setup
//#include "material_fixture.hh"

using namespace detran_test;
using namespace detran_diffusion;
using namespace detran_geometry;
using namespace detran_material;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_LossOperator_actual();

// Test of basic public interface
int test_LossOperator(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Call actual test.
  int result = test_LossOperator_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;
}


// Test of basic public interface
int test_LossOperator_actual()
{

  // Create input.
  LossOperator::SP_input inp(new InputDB());
  inp->put<string>("bc_west",  "reflect");
  inp->put<string>("bc_south", "reflect");
  inp->put<string>("bc_east",  "vacuum");
  inp->put<string>("bc_north", "vacuum");
  // Create material.
  LossOperator::SP_material mat(new Material(1, 2, false));
  mat->set_sigma_t(0, 0, 0.1890);
  mat->set_sigma_t(0, 1, 1.4633);
  mat->set_sigma_s(0, 0, 0, 0.1507);
  mat->set_sigma_s(0, 0, 1, 0.0000);
  mat->set_sigma_s(0, 1, 0, 0.0380);
  mat->set_sigma_s(0, 1, 1, 1.4536);
  mat->set_diff_coef(0, 0, 1 / ( 3*0.1890));
  mat->set_diff_coef(0, 1, 1 / ( 3*1.4633));
  mat->finalize();

  // phi0 = 1 / (0.1890-0.1507)
  // phi1 = (1 + phi0 * 0.0380) / (1.4633-1.4536)

  // Create mesh.
  vec_dbl cm(2, 0.0);
  cm[1] = 100.0;
  vec_int fm(1, 3);
  vec_int mt(1, 0);
  //Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));
  Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  //Mesh3D::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create the operator.
  LossOperator M(inp, mat, mesh);
  M.display();
  return 0;
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
