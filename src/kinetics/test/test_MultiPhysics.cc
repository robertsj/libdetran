//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_MultiPhysics.cc
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 *  @brief  Test of MultiPhysics class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_MultiPhysics)

// Detran headers
#include "utilities/TestDriver.hh"
#include "kinetics/MultiPhysics.hh"

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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// Test of basic public interface
int test_MultiPhysics(int argc, char *argv[])
{
  // Create physics container
  MultiPhysics physics(2);

  const int TEMPERATURE = 0;
  const int POWER       = 1;

  // Physics variables
  MultiPhysics::vec_dbl T(5, 1.0);
  MultiPhysics::vec_dbl P(5, 0.2);

  physics.add_variable(TEMPERATURE, T);
  physics.add_variable(POWER, P);
  for (int i = 0; i < 5; ++i)
  {
    TEST(soft_equiv(physics.variable(TEMPERATURE)[i], 1.0));
    TEST(soft_equiv(physics.variable(POWER)[i],       0.2));
    physics.variable(TEMPERATURE)[i] = 1.23;
    TEST(soft_equiv(physics.variable(TEMPERATURE)[i], 1.23));
  }
  physics.display();

  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_MultiPhysics.cc
//---------------------------------------------------------------------------//
