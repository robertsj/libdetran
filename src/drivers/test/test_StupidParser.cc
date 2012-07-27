//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_StupidParser.cc
 * \author Jeremy Roberts
 * \date   Jul 14, 2011
 * \brief  Test of StupidParser
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_StupidParser)

// Detran headers
#include "TestDriver.hh"
#include "StupidParser.hh"

// System headers
#include <iostream>
#include <cstdlib>

using namespace detran_test;
using namespace detran;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

// Test definitions.

int test_StupidParser(int argc, char *tmp_argv[])
{
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


//---------------------------------------------------------------------------//
//              end of test_StupidParser.cc
//---------------------------------------------------------------------------//
