//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_ConstantSource.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of ConstantSource
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_ConstantSource)

// Detran headers
#include "TestDriver.hh"
#include "ConstantSource.hh"
#include "Mesh2D.hh"
#include "QuadrupleRange.hh"

// Setup
/* ... */

using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_angle;
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

int test_ConstantSource(int argc, char *argv[])
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  // Create quadrature.
  Quadrature::SP_quadrature quad(new QuadrupleRange(2, 2));

  TEST(mesh);
  TEST(quad);

  ConstantSource q_e(1, mesh, 1.0, quad);

  TEST(soft_equiv(q_e.source(0, 0),    1.0));
  TEST(soft_equiv(q_e.source(0, 0, 0), 1.0*inv_four_pi));

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_ConstantSource.cc
//---------------------------------------------------------------------------//
