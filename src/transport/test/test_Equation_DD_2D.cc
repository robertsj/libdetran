//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Equation_DD_2D.cc
 *  @brief Test of Equation_DD_2D
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Equation_DD_2D_basic)

#include "TestDriver.hh"
#include "Equation_DD_2D.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "angle/LevelSymmetric.hh"

using namespace detran;
using namespace detran_utilities;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_test;
using namespace std;


int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_Equation_DD_2D_basic(int argc, char *argv[])
{
  // Get mesh, material, and quadrature

  SP_mesh mesh      = mesh_2d_fixture();
  SP_material mat   = material_fixture_1g();
  LevelSymmetric::SP_quadrature q(new detran_angle::LevelSymmetric(2, 2));

  TEST(mesh);
  TEST(mat);
  TEST(q);

  // Create
  Equation_DD_2D eq(mesh, mat, q, false);

  // Setup group, octant, and angle.
  eq.setup_group(0);
  eq.setup_octant(0);
  eq.setup_angle(0);

  // Create a phi and psi vector
  Equation_DD_2D::moments_type      phi(mesh->number_cells(), 0.0);
  Equation_DD_2D::angular_flux_type psi(mesh->number_cells(), 0.0);

  // Cell sweep source [n/cm^2-s-ster]
  Equation_DD_2D::moments_type source(mesh->number_cells(), 1.0);

  // Create incident and outgoing face fluxes.  Remember, these
  // are dumb arrays, so be careful!
  Equation_DD_2D::face_flux_type psi_in;
  Equation_DD_2D::face_flux_type psi_out;
  psi_in[0] = 1.0;
  psi_in[1] = 1.0;

  // solve
  eq.solve(0,       // i
           0,       // j
           0,       // k  not actually used
           source,  // these
           psi_in,  // are
           psi_out, // passed
           phi,     // by
           psi);    // reference

  // Check the results. FINISH.
  //TEST(psi_out[0] == 0.0);
  //TEST(psi_out[1] == 0.0);
  //TEST(phi[0]     == 0.0);
  //TEST(psi[0]     == 0.0);

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Equation_DD_2D.cc
//----------------------------------------------------------------------------//
