//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Homogenization.cc
 *  @brief Test of Homogenization
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_Homogenization)       \
        FUNC(test_HomogenizeCoarseMesh)

#include "utilities/TestDriver.hh"
#include "Homogenize.hh"
#include "geometry/Mesh1D.hh"
#include "material/Material.hh"
#include "geometry/test/mesh_fixture.hh"

using namespace detran;
using namespace detran_geometry;
using namespace detran_material;
using namespace detran_utilities;
using namespace detran_test;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

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

  // Mesh              |  mat 0  |  mat1  |
  //                   0         5       10
  vec_dbl cm(3, 0.0);
  cm[1] = 5;
  cm[2] = 10;
  vec_int fm(2, 1);
  vec_int mt(2, 0);
  mt[1] = 1;
  SP_mesh mesh = detran_geometry::Mesh1D::Create(fm, cm, mt);

  // State
  State::SP_state state(new State(input, mesh));
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    (state->phi(0))[i] = 1.0;
    (state->phi(1))[i] = 1.0;
  }

  // Spectrum (over materials)
  vec2_dbl spectrum(2, vec_dbl(2, 1.0));

  // Homogenization.  This should produce the same material.
  Homogenize H(mat);

  // Loop over state and spectrum based homogenization
  for (int type = 0; type < 2; ++type)
  {
    Material::SP_material mat2;
    if (type == 0)
      mat2 = H.homogenize(state, mesh, "MATERIAL");
    else
      mat2 = H.homogenize(spectrum, "MATERIAL", mesh, "MATERIAL");

    for (int m = 0; m < 2; ++m)
    {
      for (int g = 0; g < 2; ++g)
      {
        TEST(soft_equiv(mat->sigma_t(m, g),    mat2->sigma_t(m, g)));
        TEST(soft_equiv(mat->sigma_a(m, g),    mat2->sigma_a(m, g)));
        TEST(soft_equiv(mat->sigma_s(m, 0, g), mat2->sigma_s(m, 0, g)));
        TEST(soft_equiv(mat->sigma_s(m, 1, g), mat2->sigma_s(m, 1, g)));
      }
    }

    vec_int cg(1, 2);
    Material::SP_material mat3;
    if (type == 0)
      mat3 = H.homogenize(state, mesh, "MATERIAL", cg);
    else
      mat3 = H.homogenize(spectrum, "MATERIAL", mesh, "MATERIAL", cg);

    TEST(soft_equiv(mat3->sigma_t(0, 0),     1.5));
    TEST(soft_equiv(mat3->sigma_t(1, 0),     3.0));
    TEST(soft_equiv(mat3->sigma_s(0, 0, 0),  0.5));
    TEST(soft_equiv(mat3->sigma_s(1, 0, 0),  1.0));
    TEST(soft_equiv(mat3->sigma_a(0, 0),     1.0));
    TEST(soft_equiv(mat3->sigma_a(1, 0),     2.0));
    TEST(soft_equiv(mat3->diff_coef(0, 0),   0.25));
    TEST(soft_equiv(mat3->diff_coef(1, 0),   0.125));
  }

  // Test homogenization over partial groups
  vec_int cg(1, 1);
  vec_size_t fgroups(1, 0);
  Material::SP_material mat4;
  //   only 0th fine group
  mat4 = H.homogenize(spectrum, "MATERIAL", mesh, "MATERIAL", cg, fgroups);
  TEST(soft_equiv(mat4->sigma_t(0, 0), mat->sigma_t(0, 0)));
  TEST(soft_equiv(mat4->sigma_t(1, 0), mat->sigma_t(1, 0)));
  //   only 1st fine group
  fgroups[0] = 1;
  mat4 = H.homogenize(spectrum, "MATERIAL", mesh, "MATERIAL", cg, fgroups);
  TEST(soft_equiv(mat4->sigma_t(0, 0), mat->sigma_t(0, 1)));
  TEST(soft_equiv(mat4->sigma_t(1, 0), mat->sigma_t(1, 1)));
  //   over both, in reverse
  cg[0] = 2;
  fgroups.resize(2, 0);
  fgroups[0] = 1;
  mat4 = H.homogenize(spectrum, "MATERIAL", mesh, "MATERIAL", cg, fgroups);
  TEST(soft_equiv(mat4->sigma_t(0, 0), 1.5));
  TEST(soft_equiv(mat4->sigma_t(1, 0), 3.0));

  return 0;
}

/*
 *  This tests the use of homogenization given a coarse mesh.
 */
int test_HomogenizeCoarseMesh(int argc, char *argv[])
{
  // Input
  InputDB::SP_input input = InputDB::Create();
  input->put<int>("number_groups", 1);

  // Fine mesh material
  Material::SP_material mat = Material::Create(2, 1);
  mat->set_sigma_t(0, 0, 1.0);
  mat->set_sigma_a(0, 0, 0.5);
  mat->set_sigma_t(1, 0, 2.0);
  mat->finalize();

  // Fine mesh and coarse mesh
  vec_int fm(1, 5);
  vec_dbl cm(2, 0.0); cm[1] = 100.0;
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);
  //Mesh::SP_mesh mesh = Mesh2D::Create(fm, fm, cm, cm, mt);

  CoarseMesh::SP_coarsemesh mesher(new CoarseMesh(mesh, 3));
  mesher->get_coarse_mesh()->display();
  return 0;
  vec_int mat_map(mesh->number_cells(), 0);
  State::SP_state state(new State(input, mesh));
  for (int cell = 0; cell < mesh->number_cells(); ++cell)
  {
    mat_map[cell] = cell % 2;
    state->phi(0)[cell] = 1.0 + mesh->cell_to_i(cell) +  mesh->cell_to_j(cell);
  }
  mesh->add_mesh_map("MATERIAL", mat_map);
  Homogenize H(mat);
  Material::SP_material mat2 = H.homogenize(state, mesh, "COARSEMESH");
  mat2->display();
  TEST(mat2->number_materials() == 4);

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Homogenization.cc
//----------------------------------------------------------------------------//
