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

#include "utilities/TestDriver.hh"
#include "diffusion/OneGroupLossOperator.hh"
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

int test_OneGroupLossOperator_actual();

int test_OneGroupLossOperator(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Run actual test.
  int result = test_OneGroupLossOperator_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;

}

// Test of basic public interface
int test_OneGroupLossOperator_actual()
{

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
  for (int i = 1; i < 9; i++)
    TEST(soft_equiv(y_a[i], 0.1000000000000000));
  TEST(soft_equiv(y_a[9], 0.385714285714286));
  VecRestoreArray(y, &y_a);

  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_OneGroupLossOperator.cc
//---------------------------------------------------------------------------//
