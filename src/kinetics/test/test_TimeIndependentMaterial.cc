//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_TimeIndepndentMaterial.cc
 *  @author Rabab Elzohery
 *  @date   March 19, 2020
 *  @brief  Test of TimeIndepndentMaterial class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_TimeIndependentMaterial)

// Detran headers
#include "utilities/TestDriver.hh"
#include "kinetics/KineticsMaterial.hh"
#include "kinetics/TimeIndependentMaterial.hh"



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

int test_TimeIndependentMaterial(int argc, char *argv[])
{
  // Get the 1g KineticsParameters.
  KineticsMaterial::SP_material k_mat = KineticsMaterial::Create(2, 1, 8);

  // Reference
  double lambda[] = {0.012466675909352, 0.028291721655508, 0.042524366905518, 0.133041685328204,
                     0.292467164793226, 0.666487673615332, 1.634781086226286, 3.554600925948437};

  double beta[]   = {0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 0.00060, 0.000540, 0.000152};
  double beta_total =  0.006648;

  for (int i = 0; i < 8; ++i)
  {
     k_mat->set_lambda(i, lambda[i]);
    for (int m = 0; m < 2; ++m)
    {
      k_mat->set_beta(m, i, beta[i]);
      k_mat->set_sigma_t(m,0, 1);
      k_mat->set_sigma_s(m, 0, 0,  0.5);
    }
  }

  k_mat->set_velocity(0, 2200);
  k_mat->finalize();

  // TimeIndependentMaterial
  TimeIndependentMaterial::SP_material t_mat = TimeIndependentMaterial::Create(k_mat);
  for (int i = 0; i < 8; ++i)
  {
    TEST(soft_equiv(t_mat->lambda(i), lambda[i]));

    for (int m = 0; m < 2; ++m)
    {
      TEST(soft_equiv(t_mat->beta(m, i), beta[i]));
      TEST(soft_equiv(t_mat->sigma_s(m, 0, 0), 0.5));
      TEST(soft_equiv(t_mat->beta(m, i), beta[i]));
      TEST(soft_equiv(t_mat->beta_total(m), beta_total));
    }
  }

      TEST(soft_equiv(t_mat->velocity(0), 2200.0));
  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_TimeIndependentMaterial.cc
//---------------------------------------------------------------------------//
