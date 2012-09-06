//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_PinCell.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of PinCell class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_PinCell)

// Detran headers
#include "TestDriver.hh"
#include "PinCell.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_PinCell(int argc, char *argv[])
{
  // Get the pincell and mesh
  SP_pincell pin = pincell_fixture();
  SP_mesh mesh = pin->mesh();

  // Do tests
  TEST(mesh->number_cells()     == 49);
  TEST(mesh->number_cells_x()   == 7);
  TEST(mesh->number_cells_y()   == 7);
  TEST(mesh->number_cells_z()   == 1);
  // half_pitch-0.5*L-delta;
  double L = 1.60496875 * 0.54;
  double half_pitch = 0.5 * 1.26;
  double delta = 0.176223981031790 * 0.54;
  double mod_width = half_pitch-0.5*L-delta;
  TEST(soft_equiv(mesh->dx(0),  mod_width));
  TEST(soft_equiv(mesh->dy(0),  mod_width));
  TEST(soft_equiv(mesh->dz(0),  1.0));
  vec_int map = mesh->mesh_map("MATERIAL");
  TEST(map[0] == 0);
  TEST(map[22] == 1);
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Mesh3D.cc
//---------------------------------------------------------------------------//
