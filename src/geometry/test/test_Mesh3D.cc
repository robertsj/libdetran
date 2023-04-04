//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Mesh3D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Mesh3D class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Mesh3D.hh"
#include "mesh_fixture.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

TEST(Mesh3D, Basic)
{
  // Get the mesh
  SP_mesh mesh = mesh_3d_fixture();

  // Mesh properties:
  //   cm = [ 0.0  5.0 10.0]
  //   fm = [    5    5   ]
  //   mt = [ 0 1; 1 1; 1 1; 1 1]

  // Basic mesh things.
  EXPECT_EQ(mesh->number_cells()    , 1000);
  EXPECT_EQ(mesh->number_cells_x()  , 10);
  EXPECT_EQ(mesh->number_cells_y()  , 10);
  EXPECT_EQ(mesh->number_cells_z()  , 10);
  EXPECT_EQ(mesh->dx(0)             , 1.0);
  EXPECT_EQ(mesh->dy(0)             , 1.0);
  EXPECT_EQ(mesh->dz(0)             , 1.0);
  EXPECT_EQ(mesh->dimension()       , 3);
  int cell = 5 + 5*10 + 5*100;
  EXPECT_EQ(mesh->index(5, 5, 5)    , cell);
  EXPECT_EQ(mesh->cell_to_i(cell)   , 5);
  EXPECT_EQ(mesh->cell_to_j(cell)   , 5);
  EXPECT_EQ(mesh->cell_to_k(cell)   , 5);
  EXPECT_EQ(mesh->volume(cell)      , 1.0);

  vec_int mat_map = mesh->mesh_map("MATERIAL");
  EXPECT_EQ(mat_map[0]              , 0);
  EXPECT_EQ(mat_map[600]            , 1);
}

//---------------------------------------------------------------------------//
//              end of test_Mesh3D.cc
//---------------------------------------------------------------------------//
