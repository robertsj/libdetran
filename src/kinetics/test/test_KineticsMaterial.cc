//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_KineticsMaterial.cc
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 *  @brief  Test of KineticsParameters class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_KineticsMaterial)

// Detran headers
#include "utilities/TestDriver.hh"
#include "kinetics/KineticsMaterial.hh"

// System
#include <iostream>
#include <fstream>

using namespace detran_test;
using namespace detran_utilities;
using namespace detran;
using namespace std;
using detran_utilities::soft_equiv;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_KineticsMaterial(int argc, char *argv[])
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
    TEST(soft_equiv(k->lambda(i), lambda[i]));

    for (int m = 0; m < 2; ++m)
    {
      TEST(soft_equiv(k->beta(m, i), beta[i]));
      TEST(soft_equiv(k->beta_total(m), beta_total));
    }
  }
  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_KineticsMaterial.cc
//---------------------------------------------------------------------------//
