//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_LinearExternalSource.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of LinearExternalSource class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_LinearExternalSource)

// Detran headers
#include "utilities/TestDriver.hh"
#include "kinetics/LinearExternalSource.hh"
#include "geometry/Mesh1D.hh"
#include "external_source/ConstantSource.hh"
// System
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace detran_test;
using namespace detran;
using namespace detran_external_source;
using namespace detran_utilities;
using detran_utilities::soft_equiv;


int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_LinearExternalSource(int argc, char *argv[])
{

  // Create a mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0); cm[1] = 1.0;
  vec_int mt(1, 0);
  detran_geometry::Mesh1D::SP_mesh mesh =
    std::make_shared<detran_geometry::Mesh1D>(fm, cm, mt);

  // Define the times, [0, 1, 2]
  vec_dbl times(3, 1.0);
  times[1] = 2.0;
  times[2] = 3.0;

  // Define the external sources.  Constant
  // sources of strength [0, 1, 0].  I.e., this
  // produces a tent fuction from 1 to 3, peaked
  // at 2.
  LinearExternalSource::vec_source sources(3);
  sources[0] = ConstantSource::Create(2, mesh, 0.0);
  sources[1] = ConstantSource::Create(2, mesh, 1.0);
  sources[2] = ConstantSource::Create(2, mesh, 0.0);


  // Define the time dependent source
  LinearExternalSource Q(2, mesh, times, sources);

  double t = 0.0;
  for (int i = 0; i <= 20; ++i)
  {
    Q.set_time(t);
    std::printf(" time %8.4f %16.8f \n", t, Q.source(0, 0));
    t += 0.2;
  }

  Q.set_time(1.0);
  TEST(soft_equiv(Q.source(0, 0), 0.0));
  Q.set_time(1.5);
  TEST(soft_equiv(Q.source(0, 0), 0.5));
  Q.set_time(2.0);
  TEST(soft_equiv(Q.source(0, 0), 1.0));
  Q.set_time(2.5);
  TEST(soft_equiv(Q.source(0, 0), 0.5));
  Q.set_time(3.0);
  TEST(soft_equiv(Q.source(0, 0), 0.0));

  return 0;
}



//---------------------------------------------------------------------------//
//              end of test_KineticsMaterial.cc
//---------------------------------------------------------------------------//
