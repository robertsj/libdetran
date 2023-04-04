//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_MomentToDiscrete.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of MomentToDiscrete class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "MomentToDiscrete.hh"
#include "QuadratureFactory.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(MomentToDiscrete, Basic)
{
  QuadratureFactory qf;
  Quadrature::SP_quadrature q;
  InputDB::SP_input db = std::make_shared<InputDB>();
  // 1d test
  {
    // Indexer
    MomentIndexer::SP_momentindexer indexer = std::make_shared<MomentIndexer>(1, 0);
    // MtoD
    MomentToDiscrete MtoD(indexer);
    // Quadrature
    q = qf.build(db, 1);
    EXPECT_TRUE(q != NULL);
    // Build
    MtoD.build(q);
    // TEST
    EXPECT_EQ(MtoD(0, 0, 0), 0.5);
  }
  // 2d test
  {
    // Indexer
    MomentIndexer::SP_momentindexer indexer = std::make_shared<MomentIndexer>(2, 0);
    // MtoD
    MomentToDiscrete MtoD(indexer);
    // Quadrature
    q = qf.build(db, 2);
    EXPECT_TRUE(q != NULL);
    // Build
    MtoD.build(q);
    // TEST
    EXPECT_EQ(MtoD(0, 0, 0), inv_four_pi);
  }
}

//---------------------------------------------------------------------------//
//              end of test_MomentToDiscrete.cc
//---------------------------------------------------------------------------//
