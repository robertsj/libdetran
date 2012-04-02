//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Mesh2D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Mesh2D class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_Mesh2D_basic)

// Detran headers
#include "TestDriver.hh"
#include "Mesh2D.hh"

// Setup
#include "mesh_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace detran_utils;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Mesh2D_basic()
{
  // Get the mesh
  SP_mesh mesh = mesh_2d_fixture();

  // Basic mesh things.
  TEST(mesh->number_cells()     == 400);
  TEST(mesh->number_cells_x()   == 20);
  TEST(mesh->number_cells_y()   == 20);
  TEST(mesh->number_cells_z()   == 1);
  TEST(mesh->dx(0)              == 1.0);
  TEST(mesh->dy(0)              == 1.0);
  TEST(mesh->dz(0)              == 1.0);
  TEST(mesh->dimension()        == 2);
  TEST(mesh->index(5, 5)        == 105);

  Mesh2D::vec_int mat_map = mesh->mesh_map("MATERIAL");
  TEST(mat_map[0]               == 0);
  TEST(mat_map[10]              == 1);

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Mesh2D.cc
//---------------------------------------------------------------------------//
