//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Equation_DD_1D.cc
 *  @brief Test of Equation_DD_1D
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Equation_DD_1D)

#include "TestDriver.hh"
#include "Equation_DD_1D.hh"
#include "utilities/Definitions.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "angle/PolarQuadrature.hh"
#include <cstdio>

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

int test_Equation_DD_1D(int argc, char *argv[])
{
  // Mat
  SP_material mat(new Material(1, 1));
  mat->set_sigma_t(0, 0, 1.0);
  mat->finalize();
  // Mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0); cm[1] = 10.0;
  vec_int mt(1, 0);
  SP_mesh mesh(new Mesh1D(fm, cm, mt));
  // Quad
  PolarGL::SP_quadrature q(new PolarGL(2));
  q->display();
  // Create
  Equation_DD_1D eq(mesh, mat, q, true);

  // Setup group, octant, and angle.
  eq.setup_group(0);
  eq.setup_octant(0);
  eq.setup_angle(1);

  // Create a phi and psi vector
  Equation_DD_1D::moments_type      phi(mesh->number_cells(), 0.0);
  Equation_DD_1D::angular_flux_type psi(mesh->number_cells(), 0.0);

  // Cell sweep source [n/cm^2-s-ster]
  Equation_DD_1D::moments_type source(mesh->number_cells(), 1.0);

  // Create incident and outgoing face fluxes
  Equation_DD_1D::face_flux_type psi_in = 0.0;
  Equation_DD_1D::face_flux_type psi_out = 0.0;

  // solve
  eq.solve(0,       // i
           0,       // j  not actually used
           0,       // k  not actually used
           source,  // these
           psi_in,  // are
           psi_out, // passed
           phi,     // by
           psi);    // reference

  // Check the results. FINISH.
  printf("%20.16f %20.16f %20.16f \n", psi_out, psi[0], phi[0]);

  TEST(soft_equiv(psi_out, 0.734680275209795978));
  TEST(soft_equiv(psi[0],  0.367340137604897989));
  TEST(soft_equiv(phi[0],  0.127781046679339730));


  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Equation_DD_1D.cc
//----------------------------------------------------------------------------//
