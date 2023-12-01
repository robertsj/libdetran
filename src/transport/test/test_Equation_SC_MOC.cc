//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Equation_SC_MOC.cc
 *  @brief Test of Equation_SC_MOC
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Equation_SC_MOC.hh"
#include "Mesh2D.hh"
#include "Uniform.hh"
#include "Tracker.hh"
#include "material/test/material_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace std;

TEST(EquationSCMOC, Basic)
{

//  // Create mesh
//  vec_dbl cm(2, 0.0);
//  cm[1] = 1.0;
//  vec_int fm(1, 2);
//  vec_int mt(1, 0);
//  Mesh::SP_mesh mesh0(new Mesh2D(fm, fm, cm, cm, mt));
//
//  // Create quadrature.  1 azimuth, 3 space, 1 polar.
//  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));
//  quad->display_tracks();
//  // Create tracker
//  Tracker tracker(mesh0, quad);
//  Tracker::SP_trackdb tracks = tracker.trackdb();
//
//  // Get the tracked mesh.
//  Mesh::SP_mesh mesh = tracker.meshmoc();
//
//  // Material
//  SP_material mat = material_fixture_1g();
//
//  EXPECT_TRUE(mesh);
//  EXPECT_TRUE(mat);
//  EXPECT_TRUE(quad);
//
//  Equation_SC_MOC eq(mesh, mat, quad, true);
//
//  // Setup group, octant, and angle.
//  eq.setup_group(0);
//  eq.setup_octant(0);
//  eq.setup_azimuth(0);
//  eq.setup_polar(0);
//
//  // Create a phi and psi vector
//  Equation_SC_MOC::moments_type      phi(mesh->number_cells(), 0.0);
//  Equation_SC_MOC::angular_flux_type psi(mesh->number_cells(), 0.0);
//
//  // Cell sweep source [n/cm^2-s-ster]
//  Equation_SC_MOC::moments_type source(mesh->number_cells(), 1.0);
//
//  int region     = 0;
//  double psi_in  = 1.0;
//  double psi_out = 0.0;
//  double length  = 0.559016994;
//  double weight  = pi;
//  double space   = 0.447213595499958;
//
//  // solve
//  eq.solve(region,  // region
//           length,  // length
//           source,  // these
//           psi_in,  //   are
//           psi_out, //     passed
//           phi,     //       by
//           psi);    //         reference
//
//  // Step characteristic.
//  double t = length/0.798184;
//  double A = std::exp(-1.0*t);
//  double B = 1.0 - A;
//  double C = t * (1.0 - B/t);
//  double ref_psi_out = 1.0;
//  double ref_psi_avg = (1/t)*(B+C)*space*length/0.25;
//  double ref_phi = weight * ref_psi_avg;
//  cout << A << " " << B << " " << C << " space = " << space << endl;
//
//  cout << " psi_out " << psi_out << " ref " << ref_psi_out << endl;
//  cout << " psi_avg " << psi[0]  << " ref " << ref_psi_avg << endl;
//  cout << " psi_out " << phi[0]  << " ref " << ref_phi     << endl;
//  EXPECT_NEAR(psi_out, ref_psi_out, 1.0e-12);
//  EXPECT_NEAR(psi[0],  ref_psi_avg, 1.0e-12);
//  EXPECT_NEAR(phi[0],  ref_phi, 1.0e-12);

}

//----------------------------------------------------------------------------//
//              end of test_Equation_SC_MOC.cc
//----------------------------------------------------------------------------//
