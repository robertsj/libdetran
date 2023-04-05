//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_PinCell.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of PinCell class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "PinCell.hh"
#include "geometry/test/mesh_fixture.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

TEST(PinCell, Basic)
{
  // Get the pincell and mesh
  SP_pincell pin = pincell_fixture();
  SP_mesh mesh = pin->mesh();

  // Do tests
  EXPECT_EQ(mesh->number_cells()    , 49);
  EXPECT_EQ(mesh->number_cells_x()  , 7);
  EXPECT_EQ(mesh->number_cells_y()  , 7);
  EXPECT_EQ(mesh->number_cells_z()  , 1);
  // half_pitch-0.5*L-delta;
  double L = 1.60496875 * 0.54;
  double half_pitch = 0.5 * 1.26;
  double delta = 0.176223981031790 * 0.54;
  double mod_width = half_pitch-0.5*L-delta;
  EXPECT_NEAR(mesh->dx(0),  mod_width, 1.0e-12);
  EXPECT_NEAR(mesh->dy(0),  mod_width, 1.0e-12);
  EXPECT_NEAR(mesh->dz(0),  1.0, 1.0e-12);
  vec_int map = mesh->mesh_map("MATERIAL");
  EXPECT_EQ(map[0], 0);
  EXPECT_EQ(map[22], 1);
}

//---------------------------------------------------------------------------//
//              end of test_Mesh3D.cc
//---------------------------------------------------------------------------//
