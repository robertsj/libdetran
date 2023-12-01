//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_KineticsMaterial.cc
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 *  @brief  Test of KineticsParameters class.
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "kinetics/KineticsMaterial.hh"
#include <iostream>
#include <fstream>

using namespace detran_utilities;
using namespace detran;
using namespace std;

// Test of basic public interface
TEST(KineticsMaterial, Basic)
{
  // Get the 1g KineticsParameters.
  KineticsMaterial::SP_material k(new KineticsMaterial(2, 1, 8));

  // Reference
  double lambda[] = {0.012466675909352, 0.028291721655508, 0.042524366905518, 0.133041685328204,
                     0.292467164793226, 0.666487673615332, 1.634781086226286, 3.554600925948437};

  double beta[]   = {0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 0.00060, 0.000540, 0.000152};
  double beta_total =  0.006648;

  for (int i = 0; i < 8; ++i)
  {
    k->set_lambda(i, lambda[i]);
    for (int m = 0; m < 2; ++m)
    {
      k->set_beta(m, i, beta[i]);
    }
  }
  k->set_velocity(0, 2200);
  k->finalize();
  for (int i = 0; i < 8; ++i)
  {
    EXPECT_NEAR(k->lambda(i), lambda[i], 1.0e-12);

    for (int m = 0; m < 2; ++m)
    {
      EXPECT_NEAR(k->beta(m, i), beta[i], 1.0e-12);
      EXPECT_NEAR(k->beta_total(m), beta_total, 1.0e-12);
    }
  }
}

//---------------------------------------------------------------------------//
//              end of test_KineticsMaterial.cc
//---------------------------------------------------------------------------//
