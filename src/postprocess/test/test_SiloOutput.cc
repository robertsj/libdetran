//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_SiloOutput.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Material class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_SiloOutput)

// Detran headers
#include "TestDriver.hh"
#include "SiloOutput.hh"
#include "Mesh2D.hh"

// System
#include <iostream>

using namespace detran_test;
using namespace detran_postprocess;
using namespace detran;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

// Test of basic public interface
int test_SiloOutput(int argc, char *argv[])
{

  // Create input.
  SiloOutput::SP_input inp(new InputDB());
  inp->put<int>("number_groups", 1);
  inp->put<string>("silo_filename", "test.silo");

  // Create mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0);
  cm[1] = 10;
  vec_int mt(1, 0);
  Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  // Make state and fill.
  State::SP_state state(new State(inp, mesh));
  for (int i = 0; i < mesh->number_cells(); i++)
  {
    state->phi(0)[i] = (double) i;
  }

  // Make SiloOuput
  SiloOutput out(inp, mesh);
  out.initialize();
  out.write_flux(state);
  out.finalize();

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_SiloOutput.cc
//---------------------------------------------------------------------------//
