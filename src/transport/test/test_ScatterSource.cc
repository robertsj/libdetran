//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_ScatterSource.cc
 *  @brief Test of ScatterSource
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                        \
        FUNC(test_ScatterSource_forward) \
        FUNC(test_ScatterSource_adjoint) \

#include "utilities/TestDriver.hh"
#include "transport/ScatterSource.hh"
#include "geometry/Mesh1D.hh"

using namespace detran;
using namespace detran_geometry;
using namespace detran_material;
using namespace detran_utilities;
using namespace detran_test;
using std::cout;
using std::endl;
#define COUT(c) cout << c << endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
ScatterSource::SP_material get_material()
{
  // Material
  //
  //     |  1   2   0  |  g=0, lower=0, upper=1, lowerT=0, upperT=2
  // S = |  3   4   5  |  g=1, lower=0, upper=2, lowerT=0, upperT=2
  //     |  6   7   8  |  g=2, lower=0, upper=2, lowerT=1, upperT=2
  //
  //     |  1   3   6  |
  // S*= |  2   4   7  |
  //     |  0   5   8  |
  ScatterSource::SP_material mat(new Material(1, 3));
  double S[3][3] = {{1.0,2.0,0.0},{3.0,4.0,5.0},{6.0,7.0,8.0}};
  for (unsigned int gp = 0; gp < 3; ++gp)
    for (unsigned int g = 0; g < 3; ++g)
      mat->set_sigma_s(0, g, gp, S[g][gp]);
  mat->finalize();
  return mat;
}

//----------------------------------------------------------------------------//
State::SP_state get_state(int adjoint)
{
  // Flux is just [1, 2, 3]
  vec_int mt(1, 0);
  vec_int fm(1, 1);
  vec_dbl cm(2, 0); cm[1] = 1.0;
  ScatterSource::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);
  State::SP_input db = InputDB::Create();
  db->put<int>("number_groups", 3);
  db->put<int>("adjoint",       adjoint);
  State::SP_state state(new State(db, mesh));
  for (int g = 0; g < 3; ++g)
    state->phi(g)[0] = 1.0 + g;
  return state;
}

//----------------------------------------------------------------------------//
int test_ScatterSource_forward(int argc, char *argv[])
{
  // Setup
  ScatterSource::SP_material mat = get_material();
  State::SP_state state = get_state(0);

  // Scatter source
  ScatterSource source(state->get_mesh(), mat, state);
  vec_dbl s(1, 0.0);

  // Within-group source
  // g = 0
  source.build_within_group_source(0, state->phi(0), s);
  TEST(soft_equiv(s[0], 1.0));
  s[0] = 0.0;
  // g = 1
  source.build_within_group_source(1, state->phi(1), s);
  TEST(soft_equiv(s[0], 8.0));
  s[0] = 0.0;
  // g = 2
  source.build_within_group_source(2, state->phi(2), s);
  TEST(soft_equiv(s[0], 24.0));
  s[0] = 0.0;

  // In-scatter source
  // g = 0
  source.build_in_scatter_source(0, s);
  TEST(soft_equiv(s[0], 4.0));
  s[0] = 0.0;
  // g = 1
  source.build_in_scatter_source(1, s);
  TEST(soft_equiv(s[0], 18.0));
  s[0] = 0.0;
  // g = 2
  source.build_in_scatter_source(2, s);
  TEST(soft_equiv(s[0], 20.0));
  s[0] = 0.0;

  // Downscatter source
  // g = 0
  source.build_downscatter_source(0, 0, s);
  TEST(soft_equiv(s[0], 0.0));
  s[0] = 0.0;
  // g = 1
  source.build_downscatter_source(1, 1, s);
  TEST(soft_equiv(s[0], 3.0));
  s[0] = 0.0;
  // g = 2
  source.build_downscatter_source(2, 2, s);
  TEST(soft_equiv(s[0], 20.0));
  s[0] = 0.0;

  // Total group source
  // g = 0
  source.build_total_group_source(0, 0, state->all_phi(), s);
  TEST(soft_equiv(s[0], 5.0));
  s[0] = 0.0;
  // g = 1
  source.build_total_group_source(1, 0, state->all_phi(), s);
  TEST(soft_equiv(s[0], 26.0));
  s[0] = 0.0;
  // g = 2
  source.build_total_group_source(2, 0, state->all_phi(), s);
  TEST(soft_equiv(s[0], 44.0));
  s[0] = 0.0;

  return 0;
}

//----------------------------------------------------------------------------//
int test_ScatterSource_adjoint(int argc, char *argv[])
{
  // Setup
  ScatterSource::SP_material mat = get_material();
  State::SP_state state = get_state(1);

  // Scatter source
  ScatterSource source(state->get_mesh(), mat, state);
  vec_dbl s(1, 0.0);

  // Within-group source
  // g = 0
  source.build_within_group_source(0, state->phi(0), s);
  TEST(soft_equiv(s[0], 1.0));
  s[0] = 0.0;
  // g = 1
  source.build_within_group_source(1, state->phi(1), s);
  TEST(soft_equiv(s[0], 8.0));
  s[0] = 0.0;
  // g = 2
  source.build_within_group_source(2, state->phi(2), s);
  TEST(soft_equiv(s[0], 24.0));
  s[0] = 0.0;

  // In-scatter source
  // g = 0
  source.build_in_scatter_source(0, s);
  TEST(soft_equiv(s[0], 24.0));
  s[0] = 0.0;
  // g = 1
  source.build_in_scatter_source(1, s);
  TEST(soft_equiv(s[0], 23.0));
  s[0] = 0.0;
  // g = 2
  source.build_in_scatter_source(2, s);
  TEST(soft_equiv(s[0], 10.0));
  s[0] = 0.0;

  // Upscatter source
  // g = 0
  source.build_downscatter_source(0, 0, s);
  TEST(soft_equiv(s[0], 24.0));
  s[0] = 0.0;
  // g = 1
  source.build_downscatter_source(1, 1, s);
  TEST(soft_equiv(s[0], 21.0));
  s[0] = 0.0;
  // g = 2
  source.build_downscatter_source(2, 2, s);
  TEST(soft_equiv(s[0], 0.0));
  s[0] = 0.0;

  // Total group source
  // g = 0
  source.build_total_group_source(0, 2, state->all_phi(), s);
  TEST(soft_equiv(s[0], 25.0));
  s[0] = 0.0;
  // g = 1
  source.build_total_group_source(1, 2, state->all_phi(), s);
  TEST(soft_equiv(s[0], 31.0));
  s[0] = 0.0;
  // g = 2
  source.build_total_group_source(2, 2, state->all_phi(), s);
  TEST(soft_equiv(s[0], 34.0));
  s[0] = 0.0;

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_ScatterSource.cc
//----------------------------------------------------------------------------//
