//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Mesh1D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Mesh1D class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_Mesh1D)          \
        FUNC(test_Mesh1D_serialize)

// Detran headers
#include "TestDriver.hh"
#include "Mesh1D.hh"

// System headers
#include <fstream>

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

int test_Mesh1D(int argc, char *argv[])
{
  // Get the mesh
  SP_mesh mesh = mesh_1d_fixture();

  // Mesh properties:
  //   cm = [ 0.0  5.0 10.0]
  //   fm = [    5    5    ]
  //   mt = [    0    1    ]

  // Basic mesh things.
  TEST(mesh->number_cells()     == 10);
  TEST(mesh->number_cells_x()   == 10);
  TEST(mesh->number_cells_y()   == 1);
  TEST(mesh->number_cells_z()   == 1);
  TEST(mesh->dx(0)              == 1.0);
  TEST(mesh->dy(0)              == 1.0);
  TEST(mesh->dz(0)              == 1.0);
  TEST(mesh->dimension()        == 1);
  int cell = 5;
  TEST(mesh->index(5)           == cell);
  TEST(mesh->cell_to_i(cell)    == 5);
  TEST(mesh->cell_to_j(cell)    == 0);
  TEST(mesh->cell_to_k(cell)    == 0);
  TEST(mesh->volume(cell)       == 1.0);

  vec_int mat_map = mesh->mesh_map("MATERIAL");
  TEST(mat_map[0]               == 0);
  TEST(mat_map[5]               == 1);
  TESTFALSE(mesh->mesh_map_exists("SHOULDNOTEXIST"));

  return 0;
}

int test_Mesh1D_serialize(int argc, char *argv[])
{
#ifdef DETRAN_ENABLE_BOOST
  {
    // Get the mesh and pack it

    SP_mesh mesh = mesh_1d_fixture();
    mesh->display();
    std::ofstream ofs("mesh1d.archive");
    boost::archive::binary_oarchive oa(ofs);
    oa << mesh;
    ofs.close();
  }

  {
    // Unpack

    SP_mesh mesh;
    std::ifstream ifs("mesh1d.archive");
    boost::archive::binary_iarchive ia(ifs);
    ia >> mesh;
    ifs.close();
    mesh->display();

    // Test

    TEST(mesh->number_cells()     == 10);
    TEST(mesh->number_cells_x()   == 10);
    TEST(mesh->number_cells_y()   == 1);
    TEST(mesh->number_cells_z()   == 1);
    TEST(mesh->dx(0)              == 1.0);
    TEST(mesh->dy(0)              == 1.0);
    TEST(mesh->dz(0)              == 1.0);
    TEST(mesh->dimension()        == 1);
    int cell = 5;
    TEST(mesh->index(5)           == cell);
    TEST(mesh->cell_to_i(cell)    == 5);
    TEST(mesh->cell_to_j(cell)    == 0);
    TEST(mesh->cell_to_k(cell)    == 0);
    TEST(mesh->volume(cell)       == 1.0);
    vec_int mat_map = mesh->mesh_map("MATERIAL");
    TEST(mat_map[0]               == 0);
    TEST(mat_map[5]               == 1);
    TESTFALSE(mesh->mesh_map_exists("SHOULDNOTEXIST"));
  }
#endif
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Mesh3D.cc
//---------------------------------------------------------------------------//
