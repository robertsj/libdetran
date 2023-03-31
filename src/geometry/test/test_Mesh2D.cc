//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Mesh2D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Mesh2D class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Mesh2D.hh"
#include "mesh_fixture.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

TEST(Mesh2D, Basic)
{
  // Get the mesh
  SP_mesh mesh = mesh_2d_fixture();

  // Basic mesh things.
  EXPECT_EQ(mesh->number_cells()    , 400);
  EXPECT_EQ(mesh->number_cells_x()  , 20);
  EXPECT_EQ(mesh->number_cells_y()  , 20);
  EXPECT_EQ(mesh->number_cells_z()  , 1);
  EXPECT_EQ(mesh->dx(0)             , 1.0);
  EXPECT_EQ(mesh->dy(0)             , 1.0);
  EXPECT_EQ(mesh->dz(0)             , 1.0);
  EXPECT_EQ(mesh->dimension()       , 2);
  int cell = 5 + 5*20;
  EXPECT_EQ(mesh->index(5, 5)       , cell);
  EXPECT_EQ(mesh->cell_to_i(cell)   , 5);
  EXPECT_EQ(mesh->cell_to_j(cell)   , 5);

  vec_int mat_map = mesh->mesh_map("MATERIAL");
  EXPECT_EQ(mat_map[0]              , 0);
  EXPECT_EQ(mat_map[10]             , 1);
}

//---------------------------------------------------------------------------//
//              end of test_Mesh2D.cc
//---------------------------------------------------------------------------//
