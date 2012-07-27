//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_GainOperator.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of GainOperator class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_GainOperator)

// Detran headers
#include "TestDriver.hh"
#include "GainOperator.hh"
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

int test_GainOperator_actual();

// Test of basic public interface
int test_GainOperator(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Call actual test.
  int result = test_GainOperator_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;
}

// Test of basic public interface
int test_GainOperator_actual()
{

  // Create input.
  GainOperator::SP_input inp(new InputDB());
  inp->put<string>("bc_left",  "reflect");
  inp->put<string>("bc_right", "reflect");

  // Create material.
  GainOperator::SP_material mat(new Material(2, 1, false));
  mat->set_sigma_t(0, 0, 0.2263); // (matid, g, value);
  mat->set_sigma_t(0, 1, 1.0119);
  mat->set_sigma_f(0, 0, 0.0067);
  mat->set_sigma_f(0, 1, 0.1241);
  mat->set_chi(0, 0, 1.0);
  mat->set_chi(0, 1, 0.0);
  mat->set_sigma_s(0, 0, 0, 0.2006); // 1 <- 1
  mat->set_sigma_s(0, 0, 1, 0.0000); // 1 <- 2
  mat->set_sigma_s(0, 1, 0, 0.0161); // 2 <- 1
  mat->set_sigma_s(0, 1, 1, 0.9355); // 2 <- 2
  mat->set_diff_coef(0, 0, 0.2263/3.0);
  mat->set_diff_coef(0, 1, 1.0119/3.0);
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
  GainOperator M(inp, mat, mesh);

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

}

//---------------------------------------------------------------------------//
//              end of test_GainOperator.cc
//---------------------------------------------------------------------------//
