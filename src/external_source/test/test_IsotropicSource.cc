//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_IsotropicSource.cc
 *  @brief Test of IsotropicSource
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "IsotropicSource.hh"
#include "Mesh2D.hh"
#include "QuadratureFactory.hh"

using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(IsotropicSource, Basic)
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

  // Create spectra.  Two different spectra: the first
  // is zero (the only way to have no source with this
  // source type), and the second has a strength of
  // 1 and 2 in groups 0 and 1, respectively.
  IsotropicSource::size_t num_spectra = 2;
  IsotropicSource::size_t num_groups  = 2;
  vec2_dbl spectra(num_spectra, vec_dbl(num_groups, 0.0));
  spectra[1][0] = 1.0;
  spectra[1][1] = 2.0;
  vec_int map(mesh->number_cells(), 0);
  map[0] = 1;

  //IsotropicSource q_e(num_groups, mesh, spectra, map, quad);
  auto q_e = std::make_shared<IsotropicSource>(num_groups, mesh, spectra, map, quad);

  EXPECT_NEAR(q_e->source(0, 0),    1.0,             1.0e-12);
  EXPECT_NEAR(q_e->source(0, 0, 0), 1.0*inv_four_pi, 1.0e-12);
  EXPECT_NEAR(q_e->source(0, 1),    2.0,             1.0e-12);
  EXPECT_NEAR(q_e->source(0, 1, 0), 2.0*inv_four_pi, 1.0e-12);
  EXPECT_NEAR(q_e->source(1, 0),    0.0,             1.0e-12);
  EXPECT_NEAR(q_e->source(1, 0, 0), 0.0,             1.0e-12);
}

//----------------------------------------------------------------------------//
//              end of test_IsotropicSource.cc
//----------------------------------------------------------------------------//
