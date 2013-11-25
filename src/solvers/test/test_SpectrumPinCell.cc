//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_SpectrumPinCell.cc
 *  @brief Test of SpectrumPinCell
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                            \
        FUNC(test_SpectrumPinCell)

#include "TestDriver.hh"
#include "solvers/mg/SpectrumPinCell.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/test/fixedsource_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
int test_SpectrumPinCell(int argc, char *argv[])
{
  callow_initialize(argc, argv);

  {
    FixedSourceData data = get_fixedsource_data(1, 7);

    InputDB::SP_input db(new InputDB());
    db->put<vec_dbl>("radii_0", vec_dbl(0));
    db->put<vec_int>("mat_map_0", vec_int(1, 0));
    db->put<double>("pitch_0", 1.26);

    data.mesh->add_mesh_map("PINCELL_SPECTRUM",
                            vec_int(data.mesh->number_cells(), 0));

    data.input->put<InputDB::SP_input>("mgpc_spectrum_pincell_db", db);

    SpectrumPinCell S(data.input, data.material, data.mesh);
    SpectrumPinCell::vec2_dbl x = S.spectrum();

    double ref[] = {6.873989792058641e-02, 9.961052698673564e-01,
        5.519933833019457e-02, 1.452704956055006e-03, 2.005421939335062e-04,
        1.595459951736601e-05, 9.204673721554377e-07};

    for (int g = 0; g < 7; ++g)
    {
      TEST(soft_equiv(x[g][0], ref[g], 1.0e-6));
    }

  }
  callow_finalize();

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_SpectrumPinCell.cc
//----------------------------------------------------------------------------//
