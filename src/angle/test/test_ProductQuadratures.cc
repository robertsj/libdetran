//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_ProductQuadratures.cc
 *  @brief Test of GaussLegendre class
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_ProductQuadratures)

#include "TestDriver.hh"
#include "QuadratureFactory.hh"
#include "ProductQuadrature.hh"

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

// Make sure the a few quadratures give expected numbers
int test_ProductQuadratures(int argc, char *argv[])
{
  int N = 3;
  const char* types[] = {"u-gl", "u-dgl", "asqr-asdr"};
  double mu[3][12] = {{0.952453550758000, 0.735786494829000, 0.354307382859200,
                       0.807451131668000, 0.623769670912500, 0.300367295623300,
                       0.539521597344300, 0.416789569085700, 0.200699010465100,
                       0.189454790600400, 0.146357033572400, 0.070476120304600},
                      {0.974536571519800, 0.849384968487000, 0.452334157760700,
                       0.826172212702000, 0.720073806728800, 0.383470382661000,
                       0.552030623770000, 0.481137935381400, 0.256226717968800,
                       0.193847376539100, 0.168953174898500, 0.089974858166900},
                      {0.948499936863400, 0.635950440054300, 0.229700032293800,
                       0.827163947957100, 0.554597060322400, 0.200315865266500,
                       0.476302301137400, 0.319351268497600, 0.115347033458600,
                       0.106826722495100, 0.071625203692600, 0.025870430406300}};

  InputDB::SP_input db = InputDB::Create();
  db->put<int>("quad_number_azimuth_octant", 4);
  db->put<int>("quad_number_polar_octant",   3);
  for (int i = 0; i < N; ++i)
  {
    db->put<std::string>("quad_type", types[i]);
    QuadratureFactory::SP_quadrature q = QuadratureFactory::build(db, 3);
    for (int j = 0; j < q->number_angles_octant(); ++j)
    {
      TEST(soft_equiv(mu[i][j], q->mu(0, j), 1e-11));
    }
  }

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_ProductQuadratures.cc
//---------------------------------------------------------------------------//
