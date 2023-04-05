//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_ConstantSource.cc
 *  @brief Test of ConstantSource
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "ConstantSource.hh"
#include "Mesh2D.hh"
#include "QuadratureFactory.hh"

using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(ConstantSource, Basic)
{
  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  EXPECT_TRUE(mesh != NULL);

  // Create quadrature.
  Quadrature::SP_quadrature quad;
  QuadratureFactory qf;
  InputDB::SP_input db =std::make_shared<InputDB>();
  quad = qf.build(db, 2);
  EXPECT_TRUE(quad != NULL);

  ConstantSource q_e(1, mesh, 1.0, quad);

  EXPECT_NEAR(q_e.source(0, 0),    1.0,             1.0e-12);
  EXPECT_NEAR(q_e.source(0, 0, 0), 1.0*inv_four_pi, 1.0e-12);
}

//----------------------------------------------------------------------------//
//              end of test_ConstantSource.cc
//----------------------------------------------------------------------------//
