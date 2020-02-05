//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_SpectrumGS.cc
 *  @brief Test of SpectrumGS
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                            \
        FUNC(test_SpectrumGS)

#include "TestDriver.hh"
#include "solvers/mg/SpectrumGS.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/test/fixedsource_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
int test_SpectrumGS(int argc, char *argv[])
{
  FixedSourceData data = get_fixedsource_data(1, 7);

  {
    double ref[] = {7.58295118996886e-02, 9.95588875144603e-01,
        5.52315297625087e-02, 1.45402595100291e-03, 2.01204973835345e-04,
        1.60526070885210e-05, 9.26214237060863e-07};
    data.input->put<int>("mgpc_spectrum_include_fission", 1);
    SpectrumGS S(data.input, data.material, data.mesh);
    SpectrumGS::vec2_dbl x = S.spectrum();
    for (int g = 0; g < 7; ++g)
    {
      printf(" %20.16e  %20.16e \n", x[g][0], ref[g]);
      TEST(soft_equiv(x[g][0], ref[g], 1.0e-8));
    }
  }

  {
    double ref[] = {0, 0, 0, 1.847100932670026e-01, 9.645871496558331e-01,
        1.879791855878827e-01, 1.084610297973592e-02};
    data.input->put<int>("mgpc_spectrum_include_fission", 0);
    SpectrumGS S(data.input, data.material, data.mesh);
    SpectrumGS::vec2_dbl x = S.spectrum();
    for (int g = 0; g < 7; ++g)
    {
      printf(" %20.16e  %20.16e \n", x[g][0], ref[g]);
      TEST(soft_equiv(x[g][0], ref[g]));
    }
  }
  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_SpectrumGS.cc
//----------------------------------------------------------------------------//
