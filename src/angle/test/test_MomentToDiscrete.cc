//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_MomentToDiscrete.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of MomentToDiscrete class
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                      \
        FUNC(test_MomentToDiscrete_1D) \
        FUNC(test_MomentToDiscrete_2D)
// Detran headers
#include "TestDriver.hh"
#include "MomentToDiscrete.hh"
//
#include "GaussLegendre.hh"
#include "QuadrupleRange.hh"

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_MomentToDiscrete_1D()
{
  typedef MomentToDiscrete<_1D> M2D;

  // Get quadrature.
  M2D::SP_quadrature q;
  q = new GaussLegendre(2);

  // M2D
  M2D m(0);
  m.build(q);

  // Tests
  for (int o = 0; o < q->number_octants(); o++)
  {
    for(int a = 0; a < q->number_angles_octant(); a++)
    {
      TEST(soft_equiv(m(o, a, 0, 0), 0.5));
    }
  }

  return 0;
}

int test_MomentToDiscrete_2D()
{
  typedef MomentToDiscrete<_2D> M2D;

  // Get quadrature.
  M2D::SP_quadrature q;
  q = new QuadrupleRange(18);

  // M2D
  M2D m(0);
  m.build(q);

  // Tests
  for (int o = 0; o < q->number_octants(); o++)
  {
    for(int a = 0; a < q->number_angles_octant(); a++)
    {
      std::cout << " m2d = " << m(o, a, 0, 0) << std::endl;
      TEST(soft_equiv(m(o, a, 0, 0), inv_four_pi));
    }
  }
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_MomentToDiscrete.cc
//---------------------------------------------------------------------------//
