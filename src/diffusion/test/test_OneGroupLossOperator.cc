//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_OneGroupLossOperator.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of OneGroupLossOperator class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_OneGroupLossOperator)

// Detran headers
#include "TestDriver.hh"
#include "OneGroupLossOperator.hh"
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
int test_OneGroupLossOperator()
{

  // Initialize PETSc
  int argc = 0;
  char **args;
  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);

  // Create input.
  OneGroupLossOperator::SP_input inp(new InputDB());
  inp->put<string>("bc_left",  "vacuum");
  inp->put<string>("bc_right", "vacuum");

  // Create material.
  OneGroupLossOperator::SP_material mat(new Material(1, 1, false));
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_s(0, 0, 0, 0.9);
  mat->set_diff_coef(0, 0, 1.0/3.0);
  mat->finalize();

  // Create mesh.
  vec_dbl cm(2, 0.0);
  cm[1] = 10.0;
  vec_int fm(1, 10);
  vec_int mt(1, 0);
  Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));
  //Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  //Mesh3D::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));

  // Create the operator.
  OneGroupLossOperator M(inp, mat, mesh, 0);

  // Create a right hand side.
  Vec x;
  Vec y;
  int n = mesh->number_cells();
  VecCreateSeq(PETSC_COMM_SELF, n, &x);
  VecDuplicate(x, &y);
  VecSet(x, 1.0);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  M.multiply(x, y);
  double *y_a;
  VecGetArray(y, &y_a);
  TEST(soft_equiv(y_a[0], 0.385714285714286));
  TEST(soft_equiv(y_a[1], 0.1000000000000000));
  VecRestoreArray(y, &y_a);
 // PetscFinalize();

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_OneGroupLossOperator.cc
//---------------------------------------------------------------------------//
