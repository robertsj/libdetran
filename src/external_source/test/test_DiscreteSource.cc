//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_DiscreteSource.cc
 *  @brief Test of DiscreteSource
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "DiscreteSource.hh"
#include "Mesh2D.hh"
#include "QuadratureFactory.hh"

using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(DiscreteSource, Basic)
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

  // Create spectra.
  DiscreteSource::size_t num_spectra = 2;
  DiscreteSource::size_t num_groups  = 2;
  vec3_dbl spectra(num_spectra,
                   vec2_dbl(num_groups,
                            vec_dbl(quad->number_angles(), 0.0)));

  for (int g = 0; g < 2; g++)
    for (int o = 0; o < 4; o++)
      for (int a = 0; a < quad->number_angles_octant(); a++)
        spectra[1][g][quad->index(o, a)] = g/quad->weight(a);

  vec_int map(mesh->number_cells(), 0);
  map[0] = 1;

  DiscreteSource q_e(num_groups, mesh, spectra, map, quad);

  EXPECT_NEAR(q_e.source(0, 0),    0.0, 1.0e-12);
  EXPECT_NEAR(q_e.source(0, 0, 1), 0.0, 1.0e-12);
  EXPECT_NEAR(q_e.source(0, 1),    (double)quad->number_angles(), 1.0e-12);
  EXPECT_NEAR(q_e.source(0, 1, 1), 1/quad->weight(0),             1.0e-12);
  EXPECT_NEAR(q_e.source(1, 0),    0.0, 1.0e-12);
  EXPECT_NEAR(q_e.source(1, 0, 0), 0.0, 1.0e-12);

}

//----------------------------------------------------------------------------//
//              end of test_DiscreteSource.cc
//----------------------------------------------------------------------------//
