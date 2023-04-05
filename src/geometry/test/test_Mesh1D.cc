//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Mesh1D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Mesh1D class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Mesh1D.hh"
#include <fstream>
#include "mesh_fixture.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

TEST(Mesh1D, Basic)
{
  // Get the mesh
  SP_mesh mesh = mesh_1d_fixture();

  // Mesh properties:
  //   cm = [ 0.0  5.0 10.0]
  //   fm = [    5    5    ]
  //   mt = [    0    1    ]

  // Basic mesh things.
  EXPECT_EQ(mesh->number_cells()    , 10);
  EXPECT_EQ(mesh->number_cells_x()  , 10);
  EXPECT_EQ(mesh->number_cells_y()  , 1);
  EXPECT_EQ(mesh->number_cells_z()  , 1);
  EXPECT_EQ(mesh->dx(0)             , 1.0);
  EXPECT_EQ(mesh->dy(0)             , 1.0);
  EXPECT_EQ(mesh->dz(0)             , 1.0);
  EXPECT_EQ(mesh->dimension()       , 1);
  int cell = 5;
  EXPECT_EQ(mesh->index(5)          , cell);
  EXPECT_EQ(mesh->cell_to_i(cell)   , 5);
  EXPECT_EQ(mesh->cell_to_j(cell)   , 0);
  EXPECT_EQ(mesh->cell_to_k(cell)   , 0);
  EXPECT_EQ(mesh->volume(cell)      , 1.0);

  vec_int mat_map = mesh->mesh_map("MATERIAL");
  EXPECT_EQ(mat_map[0]              , 0);
  EXPECT_EQ(mat_map[5]              , 1);
  EXPECT_FALSE(mesh->mesh_map_exists("SHOULDNOTEXIST"));
}

TEST(Mesh1D, Serialize)
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

    EXPECT_EQ(mesh->number_cells()    , 10);
    EXPECT_EQ(mesh->number_cells_x()  , 10);
    EXPECT_EQ(mesh->number_cells_y()  , 1);
    EXPECT_EQ(mesh->number_cells_z()  , 1);
    EXPECT_EQ(mesh->dx(0)             , 1.0);
    EXPECT_EQ(mesh->dy(0)             , 1.0);
    EXPECT_EQ(mesh->dz(0)             , 1.0);
    EXPECT_EQ(mesh->dimension()       , 1);
    int cell = 5;
    EXPECT_EQ(mesh->index(5)          , cell);
    EXPECT_EQ(mesh->cell_to_i(cell)   , 5);
    EXPECT_EQ(mesh->cell_to_j(cell)   , 0);
    EXPECT_EQ(mesh->cell_to_k(cell)   , 0);
    EXPECT_EQ(mesh->volume(cell)      , 1.0);
    vec_int mat_map = mesh->mesh_map("MATERIAL");
    EXPECT_EQ(mat_map[0]              , 0);
    EXPECT_EQ(mat_map[5]              , 1);
    TESTFALSE(mesh->mesh_map_exists("SHOULDNOTEXIST"));
  }
#endif
}

//---------------------------------------------------------------------------//
//              end of test_Mesh3D.cc
//---------------------------------------------------------------------------//
