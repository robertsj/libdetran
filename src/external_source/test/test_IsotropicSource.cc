//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_IsotropicSource.cc
 *  @brief Test of IsotropicSource
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_IsotropicSource)

#include "TestDriver.hh"
#include "IsotropicSource.hh"
#include "Mesh2D.hh"
#include "QuadratureFactory.hh"

using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_IsotropicSource(int argc, char *argv[])
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
  TEST(mesh);

  // Create quadrature.
  Quadrature::SP_quadrature quad;
  QuadratureFactory qf;
  InputDB::SP_input db = InputDB::Create();
  quad = qf.build(db, 2);
  TEST(quad);

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
  IsotropicSource::SP_externalsource q_e;
  q_e = IsotropicSource::Create(num_groups, mesh, spectra, map, quad);

  TEST(soft_equiv(q_e->source(0, 0),    1.0));
  TEST(soft_equiv(q_e->source(0, 0, 0), 1.0*inv_four_pi));
  TEST(soft_equiv(q_e->source(0, 1),    2.0));
  TEST(soft_equiv(q_e->source(0, 1, 0), 2.0*inv_four_pi));
  TEST(soft_equiv(q_e->source(1, 0),    0.0));
  TEST(soft_equiv(q_e->source(1, 0, 0), 0.0));

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_IsotropicSource.cc
//----------------------------------------------------------------------------//
