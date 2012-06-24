//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Mesh3D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Mesh3D class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_Mesh3D)

// Detran headers
#include "TestDriver.hh"
#include "Mesh3D.hh"

// Setup
#include "mesh_fixture.hh"

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

int test_Mesh3D()
{
  // Get the mesh
  SP_mesh mesh = mesh_3d_fixture();

  // Mesh properties:
  //   cm = [ 0.0  5.0 10.0]
  //   fm = [    5    5   ]
  //   mt = [ 0 1; 1 1; 1 1; 1 1]

  // Basic mesh things.
  TEST(mesh->number_cells()     == 1000);
  TEST(mesh->number_cells_x()   == 10);
  TEST(mesh->number_cells_y()   == 10);
  TEST(mesh->number_cells_z()   == 10);
  TEST(mesh->dx(0)              == 1.0);
  TEST(mesh->dy(0)              == 1.0);
  TEST(mesh->dz(0)              == 1.0);
  TEST(mesh->dimension()        == 3);
  int cell = 5 + 5*10 + 5*100;
  TEST(mesh->index(5, 5, 5)     == cell);
  TEST(mesh->cell_to_i(cell)    == 5);
  TEST(mesh->cell_to_j(cell)    == 5);
  TEST(mesh->cell_to_k(cell)    == 5);
  TEST(mesh->volume(cell)       == 1.0);

  vec_int mat_map = mesh->mesh_map("MATERIAL");
  TEST(mat_map[0]               == 0);
  TEST(mat_map[600]             == 1);

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Mesh3D.cc
//---------------------------------------------------------------------------//
