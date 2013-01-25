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
#define TEST_LIST                   \
        FUNC(test_MomentToDiscrete)

// Detran headers
#include "TestDriver.hh"
#include "MomentToDiscrete.hh"
#include "QuadratureFactory.hh"

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

int test_MomentToDiscrete(int argc, char *argv[])
{
  QuadratureFactory qf;
  Quadrature::SP_quadrature q;
  InputDB::SP_input db = InputDB::Create();
  // 1d test
  {
    // Indexer
    MomentIndexer::SP_momentindexer indexer = MomentIndexer::Create(1, 0);
    // MtoD
    MomentToDiscrete MtoD(indexer);
    // Quadrature
    q = qf.build(db, 1);
    TEST(q);
    // Build
    MtoD.build(q);
    // TEST
    TEST(MtoD(0, 0, 0) == 0.5);
  }
  // 2d test
  {
    // Indexer
    MomentIndexer::SP_momentindexer indexer = MomentIndexer::Create(2, 0);
    // MtoD
    MomentToDiscrete MtoD(indexer);
    // Quadrature
    q = qf.build(db, 2);
    TEST(q);
    // Build
    MtoD.build(q);
    // TEST
    TEST(MtoD(0, 0, 0) == inv_four_pi);
  }
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_MomentToDiscrete.cc
//---------------------------------------------------------------------------//
