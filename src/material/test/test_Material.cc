//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Material.cc
 *  @brief Test of Material class.
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Material_basic)     \
        FUNC(test_Material_bounds)    \
        FUNC(test_Material_serialize)

// Detran headers
#include "TestDriver.hh"
#include "Material.hh"
#include "SoftEquivalence.hh"

// System
#include <iostream>
#include <fstream>

// Setup
#include "material_fixture.hh"

using namespace detran_test;
using namespace detran_material;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
int test_Material_basic(int argc, char *argv[])
{
  // Get the 1g material.
  SP_material mat_1g = material_fixture_1g();

  // Make sure something is there.
  TEST(mat_1g);
  TEST(mat_1g->sigma_t(0, 0)        == 1.0);
  TEST(mat_1g->sigma_s(0, 0, 0)     == 0.9);

  TEST(mat_1g->nu_sigma_f(0, 0)     == 0.0);
  TEST(mat_1g->chi(0, 0)            == 0.0);
  // add absorption ??
  TEST(mat_1g->number_groups()      == 1);
  TEST(mat_1g->number_materials()   == 3);

  // Get the 2g material
  SP_material mat_2g = material_fixture_2g();
  TEST(mat_2g);
  TEST(mat_2g->sigma_t(1, 0)        == 0.2263);
  TEST(mat_2g->sigma_s(1, 0, 0)     == 0.2006);
  TEST(mat_2g->sigma_s(1, 0, 1)     == 0.0000);
  TEST(mat_2g->sigma_s(1, 1, 0)     == 0.0161);
  TEST(mat_2g->nu_sigma_f(1, 0)     == 0.0067);
  TEST(mat_2g->chi(1, 0)            == 1.0);
  TEST(mat_2g->number_groups()      == 2);
  TEST(mat_2g->number_materials()   == 4);
  //TEST(mat_2g->downscatter()        == false);

  // Get the 7g material
  SP_material mat_7g = material_fixture_7g();
  // Add up chi for material 0; ensure = 1.0.
  double chi = 0;
  for (int g = 0; g < 7; g++)
  {
    chi += mat_7g->chi(0, g);
  }
  // false for now (0.999918)
  //TEST(soft_equiv(chi, 1.0)         == true);

  return 0;
}

// Test of group bounds for upscatter
int test_Material_bounds(int argc, char *argv[])
{

  // Get the 7g material.
  SP_material mat = material_fixture_7g();

  // Test scattering bounds
  TEST(mat->lower(0)       == 0);
  TEST(mat->upper(0)       == 0);
  TEST(mat->lower(0, true) == 4);
  TEST(mat->upper(0, true) == 0);

  TEST(mat->lower(1) == 0);
  TEST(mat->upper(1) == 1);
  TEST(mat->lower(1, true) == 6);
  TEST(mat->upper(1, true) == 1);

  TEST(mat->lower(2) == 0);
  TEST(mat->upper(2) == 2);
  TEST(mat->lower(2, true) == 6);
  TEST(mat->upper(2, true) == 2);

  TEST(mat->lower(3) == 0);
  TEST(mat->upper(3) == 4);
  TEST(mat->lower(3, true) == 6);
  TEST(mat->upper(3, true) == 3);

  TEST(mat->lower(4) == 0);
  TEST(mat->upper(4) == 5);
  TEST(mat->lower(4, true) == 6);
  TEST(mat->upper(4, true) == 3);

  TEST(mat->lower(5) == 1);
  TEST(mat->upper(5) == 6);
  TEST(mat->lower(5, true) == 6);
  TEST(mat->upper(5, true) == 4);

  TEST(mat->lower(6) == 1);
  TEST(mat->upper(6) == 6);
  TEST(mat->lower(6, true) == 6);
  TEST(mat->upper(6, true) == 5);

  // Test cutoff;
  TEST(mat->upscatter_cutoff(false) == 3);
  TEST(mat->upscatter_cutoff(true)  == 6);

  return 0;
}

// Test of basic public interface
int test_Material_serialize(int argc, char *argv[])
{
#ifdef DETRAN_ENABLE_BOOST
  // Create and pack
  {
    SP_material mat_1g = material_fixture_1g();
    std::ofstream ofs("material.archive");
    boost::archive::binary_oarchive oa(ofs);
    oa << mat_1g;
    ofs.close();
  }

  // Unpack and test
  {

    SP_material mat_1g;
    std::ifstream ifs("material.archive");
    boost::archive::binary_iarchive ia(ifs);
    ia >> mat_1g;
    ifs.close();

    TEST(mat_1g);
    TEST(soft_equiv(mat_1g->sigma_t(0, 0),    1.0));
    TEST(soft_equiv(mat_1g->sigma_s(0, 0, 0), 0.9));
    TEST(soft_equiv(mat_1g->nu_sigma_f(0, 0), 0.0));
    TEST(soft_equiv(mat_1g->chi(0, 0),        0.0));
    TEST(mat_1g->number_groups()           == 1);
    TEST(mat_1g->number_materials()        == 3);

  }
#endif
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_Material.cc
//---------------------------------------------------------------------------//
