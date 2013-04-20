//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ReactionRates.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Material class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_ReactionRates_pinpower)

// Detran headers
#include "TestDriver.hh"
#include "ReactionRates.hh"

// Fixtures
#include "material/test/material_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
//#include "solvers/test/eigenproblem_fixture.hh"

// System
#include <iostream>

using namespace detran_test;
using namespace detran_utilities;
using namespace detran_geometry;
using namespace detran_material;
using namespace detran_postprocess;
using namespace detran;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//


// Tests computation of pin power for a very simple mesh.
int test_ReactionRates_pinpower(int argc, char *argv[])
{
  // Materials
  Material::SP_material mat(new Material(1, 1));
  mat->set_sigma_f(0, 0, 1.0);
  mat->finalize();

  // Mesh
  vec_dbl cm(4, 0.0);
  cm[1] = 1.0;
  cm[2] = 2.0;
  cm[3] = 3.0;
  vec_int fm(3, 1);
  vec_int mt(9, 0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  vec_int pin_map(9, 0);
  for (int p = 0; p < 9; p++)
    pin_map[p] = p;
  mesh->add_coarse_mesh_map("PINS", pin_map);

  // Input
  InputDB::SP_input inp(new InputDB());
  inp->put<int>("number_groups", 2);

  // State
  State::SP_state state(new State(inp, mesh));
  for (int g = 0; g < 2; g++)
  {
    for (int cell = 0; cell < mesh->number_cells(); cell++)
    {
      state->phi(g)[cell] = (double) pin_map[cell];
    }
  }

  // Reaction rate
  ReactionRates rates(mat, mesh, state);

  // Compute the rates.
  vec_dbl pinpower = rates.region_power("PINS", 36.0);

  // Test
  TEST(pinpower.size() == 9);
  double sum_pinpower = 0;
  for (int i = 0; i < 9; i++)
  {
    TEST(soft_equiv(pinpower[i], (double) i));
    sum_pinpower += pinpower[i];
  }
  TEST(soft_equiv(sum_pinpower, 36.0));


  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_ReactionRates.cc
//---------------------------------------------------------------------------//
