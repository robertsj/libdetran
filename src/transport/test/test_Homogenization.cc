//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_Homogenization.cc
 *  @author Jeremy Roberts
 *  @date   Mar 27, 2013
 *  @brief  Test of Homogenization.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Homogenization)

#include "utilities/TestDriver.hh"
#include "Homogenize.hh"
#include "geometry/Mesh1D.hh"
#include "material/Material.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"
#include "angle/test/quadrature_fixture.hh"

using namespace detran;
using namespace detran_material;
using namespace detran_utilities;
using namespace detran_test;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Homogenization(int argc, char *argv[])
{

  // Input
  InputDB::SP_input input = InputDB::Create();
  input->put<int>("number_groups",    2);

  // Material
  Material::SP_material mat = Material::Create(2, 2);
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_t(0, 1, 2.0);
  mat->set_sigma_s(0, 0, 0, 0.1);
  mat->set_sigma_s(0, 1, 0, 0.2);
  mat->set_sigma_s(0, 0, 1, 0.3);
  mat->set_sigma_s(0, 1, 1, 0.4);
  mat->set_sigma_t(1, 0, 2.0);
  mat->set_sigma_t(1, 1, 4.0);
  mat->set_sigma_s(1, 0, 0, 0.2);
  mat->set_sigma_s(1, 1, 0, 0.4);
  mat->set_sigma_s(1, 0, 1, 0.6);
  mat->set_sigma_s(1, 1, 1, 0.8);
  mat->compute_sigma_a();
  mat->compute_diff_coef();
  mat->finalize();
//  mat->display();

  // Mesh              |  mat 0  |  mat1  |
  //                   0         5       10
  vec_dbl cm(3, 0.0);
  cm[1] = 5;
  cm[2] = 10;
  vec_int fm(2, 1);
  vec_int mt(2, 0);
  mt[1] = 1;
  SP_mesh mesh = detran_geometry::Mesh1D::Create(fm, cm, mt);
//  mesh->display();

  // State
  State::SP_state state(new State(input, mesh));
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    (state->phi(0))[i] = 1.0;
    (state->phi(1))[i] = 1.0;
  }

  // Homogenization.  This should produce the same material.
  Homogenize H(mat);
  Material::SP_material mat2 = H.homogenize(state, mesh, "MATERIAL");

  for (int m = 0; m < 2; ++m)
  {
    for (int g = 0; g < 2; ++g)
    {
      TEST(soft_equiv(mat->sigma_t(m, g), mat2->sigma_t(m, g)));
      TEST(soft_equiv(mat->sigma_a(m, g), mat2->sigma_a(m, g)));
      TEST(soft_equiv(mat->sigma_s(m, 0, g), mat2->sigma_s(m, 0, g)));
      TEST(soft_equiv(mat->sigma_s(m, 1, g), mat2->sigma_s(m, 1, g)));
    }
  }

  vec_int cg(1, 2);
  Material::SP_material mat3 = H.homogenize(state, mesh, "MATERIAL", cg);
//  mat3->display();

  TEST(soft_equiv(mat3->sigma_t(0, 0),     1.5));
  TEST(soft_equiv(mat3->sigma_t(1, 0),     3.0));
  TEST(soft_equiv(mat3->sigma_s(0, 0, 0),  0.5));
  TEST(soft_equiv(mat3->sigma_s(1, 0, 0),  1.0));
  TEST(soft_equiv(mat3->sigma_a(0, 0),     1.0));
  TEST(soft_equiv(mat3->sigma_a(1, 0),     2.0));
  TEST(soft_equiv(mat3->diff_coef(0, 0),   0.25));
  TEST(soft_equiv(mat3->diff_coef(1, 0),   0.125));

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_State.cc
//---------------------------------------------------------------------------//
