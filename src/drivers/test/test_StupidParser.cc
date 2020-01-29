//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_StupidParser.cc
 *  @brief  Test of StupidParser
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_StupidParser)       \
        FUNC(test_StupidParser_hdf5)

#include "utilities/TestDriver.hh"
#include "drivers/StupidParser.hh"
#include <iostream>
#include <cstdlib>

using namespace detran_test;
using namespace detran;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

static void write_test_input()
{
  using std::cout;
  using std::endl;
  fstream f;
  f.open("test.inp");
  f <<
  "# examples/simple_box.inp " << endl <<
  "# " << endl <<
  "# test input " << endl <<
  " " << endl <<
  "# Input " << endl <<
  "int number_groups   2 " << endl <<
  "int dimension       2 " << endl <<
  "str problem_type    eigenvalue " << endl <<
  "str equation        sc " << endl <<
  "int inner_max_iters 1000 " << endl <<
  "dbl inner_tolerance 1e-6 " << endl <<
  "int outer_max_iters 1000 " << endl <<
  "dbl outer_tolerance 1e-6 " << endl <<
  "int eigen_max_iters 10 " << endl <<
  "dbl eigen_tolerance 1e-4 " << endl <<
  "str bc_left         reflect " << endl <<
  "str bc_right        reflect " << endl <<
  "str bc_bottom       reflect " << endl <<
  "str bc_top          reflect " << endl <<
  "vec_int testvecint  1 2 3 4 5 " << endl <<
  "str quad_type       quadruplerange " << endl <<
  "int quad_order      50 " << endl <<
  " " << endl <<
  "# Mesh " << endl <<
  "vec_dbl mesh_xcme    0.0 2.0 4.0 " << endl <<
  "vec_dbl mesh_ycme    0.0 2.0 4.0 " << endl <<
  "vec_int mesh_xfm     10 10 " << endl <<
  "vec_int mesh_yfm     20 20 " << endl <<
  "vec_int mesh_mat_map 1 1 1 1 " << endl <<
  " " << endl <<
  "# Material " << endl <<
  "material number_materials 2 " << endl <<
  "# moderator " << endl <<
  "material sigma_t    0   0.1890 1.4633 " << endl <<
  "material sigma_s    0 0 0.1507 0.0000 " << endl <<
  "material sigma_s    0 1 0.0380 1.4536 " << endl <<
  "# fuel " << endl <<
  "material sigma_t    1   0.2263 1.0119 " << endl <<
  "material sigma_s    1 0 0.2006 0.0000 " << endl <<
  "material sigma_s    1 1 0.0161 0.9355 " << endl <<
  "material sigma_f    1   0.0067 0.1241 " << endl <<
  "material chi        1   1.0000 0.0000 " << endl <<
  " " << endl <<
  "# Fixed source " << endl <<
  "str source_type             isotropic " << endl <<
  "int source_number_spectra   2 " << endl <<
  "vec_dbl source_spectra      1.00 1.00   0.00 0.00 " << endl <<
  "vec_int source_spectra_map  0 0 0 0 " << endl;
  f.close();
}

int test_StupidParser(int argc, char *tmp_argv[])
{
  write_test_input();
  cout << endl;
  char *argv[] = {"n/a", "test.inp"};
  StupidParser parser(2, argv);

  // Parse the input.
  StupidParser::SP_input input = parser.parse_input();
  TEST(input);
  input->display();
  TEST(input->get<int>("number_groups") == 2);
  TEST(input->get<string>("bc_bottom") == "reflect");

  // Parse the material.
  StupidParser::SP_material mat = parser.parse_material();
  mat->display();
  TEST(mat->number_materials() == 2);
  TEST(mat->number_groups() == 2);
  TEST(soft_equiv(mat->sigma_t(0, 0),    0.1890));
  TEST(soft_equiv(mat->sigma_t(0, 1),    1.4633));
  TEST(soft_equiv(mat->sigma_s(0, 0, 0), 0.1507));
  TEST(soft_equiv(mat->sigma_s(0, 0, 1), 0.0000));
  TEST(soft_equiv(mat->sigma_s(0, 1, 0), 0.0380));
  TEST(soft_equiv(mat->sigma_s(0, 1, 1), 1.4536));
  TEST(soft_equiv(mat->nu_sigma_f(0, 0), 0.0000));
  TEST(soft_equiv(mat->nu_sigma_f(0, 1), 0.0000));
  TEST(soft_equiv(mat->chi(0, 0),        0.0000));
  TEST(soft_equiv(mat->chi(0, 1),        0.0000));

  TEST(soft_equiv(mat->sigma_t(1, 0),    0.2263));
  TEST(soft_equiv(mat->sigma_t(1, 1),    1.0119));
  TEST(soft_equiv(mat->sigma_s(1, 0, 0), 0.2006));
  TEST(soft_equiv(mat->sigma_s(1, 0, 1), 0.0000));
  TEST(soft_equiv(mat->sigma_s(1, 1, 0), 0.0161));
  TEST(soft_equiv(mat->sigma_s(1, 1, 1), 0.9355));
  TEST(soft_equiv(mat->sigma_f(1, 0),    0.0067));
  TEST(soft_equiv(mat->sigma_f(1, 1),    0.1241));
  TEST(soft_equiv(mat->nu(1, 0),         1.0000));
  TEST(soft_equiv(mat->nu(1, 1),         1.0000));
  TEST(soft_equiv(mat->nu_sigma_f(1, 0), 0.0067));
  TEST(soft_equiv(mat->nu_sigma_f(1, 1), 0.1241));
  TEST(soft_equiv(mat->chi(1, 0),        1.0000));
  TEST(soft_equiv(mat->chi(1, 1),        0.0000));
  // Parse the mesh.
  StupidParser::SP_mesh mesh = parser.parse_mesh();
  TEST(mesh->number_cells() == 800);

  // Check that when we insert, it shows up, and we can get it.
  return 0;
}

int test_StupidParser_hdf5(int argc, char *tmp_argv[])
{
  cout << endl;
  char *argv[] = {"n/a", "test.h5"};
  StupidParser parser(2, argv);
  StupidParser::SP_input input = parser.parse_input();
  input->display();
  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_StupidParser.cc
//---------------------------------------------------------------------------//
