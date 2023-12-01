//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_MultiPhysics.cc
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 *  @brief  Test of MultiPhysics class.
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "kinetics/MultiPhysics.hh"
#include <iostream>
#include <fstream>

using namespace detran_utilities;
using namespace detran;
using namespace std;

// Test of basic public interface
TEST(MultiPhysics, Basic)
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
    EXPECT_NEAR(physics.variable(TEMPERATURE)[i], 1.0, 1.0e-12);
    EXPECT_NEAR(physics.variable(POWER)[i],       0.2, 1.0e-12);
    physics.variable(TEMPERATURE)[i] = 1.23;
    EXPECT_NEAR(physics.variable(TEMPERATURE)[i], 1.23, 1.0e-12);
  }
  physics.display();
}

//---------------------------------------------------------------------------//
//              end of test_MultiPhysics.cc
//---------------------------------------------------------------------------//
