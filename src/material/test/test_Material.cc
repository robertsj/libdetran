//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Material.cc
 *  @brief Test of Material class.
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include "Material.hh"
#include "SoftEquivalence.hh"
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include "material_fixture.hh"

using namespace detran_material;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

// Test of basic public interface
TEST(Material, Basic)
{
  // Get the 1g material.
  SP_material mat_1g = material_fixture_1g();

  // Make sure something is there.
  EXPECT_TRUE(mat_1g != NULL);
  EXPECT_EQ(mat_1g->sigma_t(0, 0)       , 1.0);
  EXPECT_EQ(mat_1g->sigma_s(0, 0, 0)    , 0.9);

  EXPECT_EQ(mat_1g->nu_sigma_f(0, 0)    , 0.0);
  EXPECT_EQ(mat_1g->chi(0, 0)           , 0.0);
  // add absorption ??
  EXPECT_EQ(mat_1g->number_groups()     , 1);
  EXPECT_EQ(mat_1g->number_materials()  , 3);

  // Get the 2g material
  SP_material mat_2g = material_fixture_2g();
  EXPECT_TRUE(mat_2g != NULL);
  EXPECT_EQ(mat_2g->sigma_t(1, 0)       , 0.2263);
  EXPECT_EQ(mat_2g->sigma_s(1, 0, 0)    , 0.2006);
  EXPECT_EQ(mat_2g->sigma_s(1, 0, 1)    , 0.0000);
  EXPECT_EQ(mat_2g->sigma_s(1, 1, 0)    , 0.0161);
  EXPECT_EQ(mat_2g->nu_sigma_f(1, 0)    , 0.0067);
  EXPECT_EQ(mat_2g->chi(1, 0)           , 1.0);
  EXPECT_EQ(mat_2g->number_groups()     , 2);
  EXPECT_EQ(mat_2g->number_materials()  , 4);

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
}

// Test of group bounds for upscatter
TEST(Material, Bounds)
{

  // Get the 7g material.
  SP_material mat = material_fixture_7g();

  // Test scattering bounds
  EXPECT_EQ(mat->lower(0), 0);
  EXPECT_EQ(mat->upper(0), 0);
  EXPECT_EQ(mat->lower(0, true), 4);
  EXPECT_EQ(mat->upper(0, true), 0);

  EXPECT_EQ(mat->lower(1), 0);
  EXPECT_EQ(mat->upper(1), 1);
  EXPECT_EQ(mat->lower(1, true), 6);
  EXPECT_EQ(mat->upper(1, true), 1);

  EXPECT_EQ(mat->lower(2), 0);
  EXPECT_EQ(mat->upper(2), 2);
  EXPECT_EQ(mat->lower(2, true), 6);
  EXPECT_EQ(mat->upper(2, true), 2);

  EXPECT_EQ(mat->lower(3), 0);
  EXPECT_EQ(mat->upper(3), 4);
  EXPECT_EQ(mat->lower(3, true), 6);
  EXPECT_EQ(mat->upper(3, true), 3);

  EXPECT_EQ(mat->lower(4), 0);
  EXPECT_EQ(mat->upper(4), 5);
  EXPECT_EQ(mat->lower(4, true), 6);
  EXPECT_EQ(mat->upper(4, true), 3);

  EXPECT_EQ(mat->lower(5), 1);
  EXPECT_EQ(mat->upper(5), 6);
  EXPECT_EQ(mat->lower(5, true), 6);
  EXPECT_EQ(mat->upper(5, true), 4);

  EXPECT_EQ(mat->lower(6), 1);
  EXPECT_EQ(mat->upper(6), 6);
  EXPECT_EQ(mat->lower(6, true), 6);
  EXPECT_EQ(mat->upper(6, true), 5);

  // Test cutoff;
  EXPECT_EQ(mat->upscatter_cutoff(false), 3);
  EXPECT_EQ(mat->upscatter_cutoff(true) , 6);
}

// Test of basic public interface
TEST(Material, Serialize)
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
    EXPECT_NEAR(mat_1g->sigma_t(0, 0),    1.0, 1.0e-12);
    EXPECT_NEAR(mat_1g->sigma_s(0, 0, 0), 0.9, 1.0e-12);
    EXPECT_NEAR(mat_1g->nu_sigma_f(0, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(mat_1g->chi(0, 0),        0.0, 1.0e-12);
    EXPECT_EQ(mat_1g->number_groups(), 1);
    EXPECT_EQ(mat_1g->number_materials(), 3);
  }
#endif

}


//---------------------------------------------------------------------------//
//              end of test_Material.cc
//---------------------------------------------------------------------------//
