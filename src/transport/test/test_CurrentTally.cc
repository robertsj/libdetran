//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_CurrentTally.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Equation_DD_2D
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_CurrentTally_1D)

// Detran headers
#include "TestDriver.hh"
#include "CurrentTally.hh"
#include "GaussLegendre.hh"
#include "Definitions.hh"

// Setup
#include "coarsemesh_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_CurrentTally_1D(int argc, char *argv[])
{
  using detran::u_int;
  typedef CurrentTally<_1D> CurrentTally_T;
  /*
   *  The 1D sample problem has a coarse mesh like
   *    0         1      2   3   4   5   6   7
   *    |    3    |   2  | 1 | 1 | 1 | 1 | 1 |
   *
   *  1001010101010101
   *
   *  where the numbers are in cm.  This was based on a fine mesh with two
   *  5 cm regions, the first of which had 5 fine meshes (now 2 of unequal
   *  widths), and the second which had 10 fine meshes (for 5 equal width
   *  coarse meshes).  Hence, the "level" factor was 2.
   *
   *  In this 1D example, we will use an S4 quadrature, so two angles
   *  per direction.  We will assume that all angular fluxes are
   *  attenuated by a factor dx, where dx is the fine mesh dx.
   *
   * We define the incident flux to be 100*o + 10*a + 1, or [1, 11, 101, and 111]
   * Using the quadrature parameters, we find that the partial currents
   * at the coarse mesh boundaries are
   *         RIGHT                  LEFT
   * 0 (2.73843733182052e+00) 5.35792775340790e-02
   * 1  2.73843733182052e+00  5.35792775340790e-02
   * 2  2.73843733182052e+00  5.35792775340790e-02
   * 3  6.84609332955131e-01  2.14317110136316e-01
   * 4  1.71152333238783e-01  8.57268440545264e-01
   * 5  4.27880833096957e-02  3.42907376218105e+00
   * 6  1.06970208274239e-02  1.37162950487242e+01
   * 7  2.67425520685598e-03 (5.48651801948969e+01)
   *
   *
   */

  // Reference partial currents
  double Jright[] = {2.73843733182052E+00, 2.73843733182052E+00, 2.73843733182052E+00,
                     6.84609332955131E-01, 1.71152333238783E-01, 4.27880833096957E-02,
                     1.06970208274239E-02, 2.67425520685598E-03};
  double Jleft[]  = {5.35792775340790E-02, 5.35792775340790E-02, 5.35792775340790E-02,
                     2.14317110136316E-01, 8.57268440545264E-01, 3.42907376218105E+00,
                     1.37162950487242E+01, 5.48651801948969E+01};

  // Get the coarse mesh
  CurrentTally_T::SP_coarsemesh mesh = coarsemesh_1d();
  TEST(mesh);

  // Get the fine mesh.
  CoarseMesh::SP_mesh finemesh = mesh->get_fine_mesh();
  TEST(finemesh->dimension() == 1);
  TEST(finemesh->number_cells() == 15);

  // Create an S4 quadrature.
  CurrentTally_T::SP_quadrature quad(new GaussLegendre(4));

  // Create the tally.
  CurrentTally_T::SP_currenttally tally(new CurrentTally_T(mesh, quad, 1));

  // Create the face flux.  For 1D, it's just a double.
  CurrentTally_T::face_flux_type psi_in, psi_out;

  // Now, fake a sweep.

  // Loop through all octants.
  for (u_int o = 0; o < 2; o++)
  {

    // Loop through azimuths in an octant.
    for (u_int a = 0; a < 2; a++)
    {

      // Start with a fixed psi at the boundary.
      psi_out = 100.0 * (double) o + 10.0 * (double) a + 1.0;

      // TALLY THE INCIDENT BOUNDARY
      u_int io = 0;
      u_int zero = 0;
      if (o == 1) io = 14;

      tally->tally(io, zero, zero, zero, o, a, psi_out, true);

      // Loop over cells.
      for (u_int ii = 0; ii < 15; ii++)
      {
        u_int i = ii;
        if (o == 1) i = 15 - i - 1;

        psi_in = psi_out;
        psi_out = psi_in * finemesh->dx(i);

        // TALLY THE OUTGOING CELL FLUX
        tally->tally(i, 0, 0, 0, o, a, psi_out, false);
      }
    }
  }

  // Test
  for (int i = 0; i < 8; i++)
  {
    TEST(soft_equiv(Jright[i], tally->partial_current(i, 0, 0, 0, 0, true)));
    TEST(soft_equiv(Jleft[i],  tally->partial_current(i, 0, 0, 0, 0, false)));
  }

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_CoarseMesh.cc
//---------------------------------------------------------------------------//
