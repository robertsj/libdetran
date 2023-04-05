//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_PulsedExternalSource.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of PulsedExternalSource class.
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "kinetics/PulsedExternalSource.hh"
#include "geometry/Mesh1D.hh"
#include "external_source/ConstantSource.hh"
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace detran;
using namespace detran_external_source;
using namespace detran_utilities;

// Test of basic public interface
TEST(PulsedExternalSource, Basic)
{

  // Create a mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0); cm[1] = 1.0;
  vec_int mt(1, 0);
  detran_geometry::Mesh1D::SP_mesh mesh =
    std::make_shared<detran_geometry::Mesh1D>(fm, cm, mt);

  // Define the fixed external source
  PulsedExternalSource::SP_externalsource fixed_source;
  fixed_source = std::make_shared<ConstantSource>(2, mesh, 1.0);

  // Define the time dependent source.  It's a pulse centered
  // at t = 2.
  PulsedExternalSource Q(2, mesh, fixed_source, 2.0, 0.2);

  double t = 0.0;
  for (int i = 0; i <= 20; ++i)
  {
    Q.set_time(t);
    std::printf(" time %8.4f %16.8f \n", t, Q.source(0, 0));
    t += 0.2;
  }


  Q.set_time(1.8);
  EXPECT_NEAR(Q.source(0, 0), 0.0625, 1.0e-12);
  Q.set_time(2.0);
  EXPECT_NEAR(Q.source(0, 0), 1.0, 1.0e-12);
  Q.set_time(2.2);
  EXPECT_NEAR(Q.source(0, 0), 0.0625, 1.0e-12);
}

//---------------------------------------------------------------------------//
//              end of test_KineticsMaterial.cc
//---------------------------------------------------------------------------//
