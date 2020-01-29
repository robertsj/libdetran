//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  test_SiloOutput.cc
 *  @brief Test of SiloOutput
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_SiloOutput)           \
        FUNC(test_SiloOutput_directory)

// Detran headers
#include "utilities/TestDriver.hh"
#include "ioutils/SiloOutput.hh"
#include "geometry/Mesh2D.hh"
#include "angle/QuadratureFactory.hh"
#include <iostream>

using namespace detran_test;
using namespace detran_utilities;
using namespace detran_geometry;
using namespace detran_angle;
using namespace detran_material;
using namespace detran_ioutils;
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
  inp->put<int>("number_groups",        1);
  inp->put<int>("store_angular_flux",   1);

  // Create mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0);
  cm[1] = 10;
  vec_int mt(1, 0);
  Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  // Quadrature
  QuadratureFactory qf;
  inp->put<int>("quad_number_polar_octant", 2);
  Quadrature::SP_quadrature quad = qf.build(inp, 2);

  // Make state and fill.
  State::SP_state state(new State(inp, mesh, quad));
  for (int i = 0; i < mesh->number_cells(); i++)
  {
    state->phi(0)[i] = (double) i;
    for (int o = 0; o < 4; o++)
      for (int a = 0; a < 2; a++)
        state->psi(0, o, a)[i] =  1000.0 * o + 100.0 * a + 1.0 * i;
  }

  // Create the SiloOutput.
  SiloOutput out(mesh);
  cout << "...created out." << endl;

  // Open the file, and write the mesh.
  TEST(out.initialize("test.silo"));
  cout << "...initialized out." << endl;

  // Write out a mesh map.
  TEST(out.write_mesh_map("MATERIAL"));
  cout << "...wrote mesh map." << endl;

  // Write out the scalar flux.
  TEST(out.write_scalar_flux(state));
  cout << "...wrote scalar flux." << endl;

  // Write out the angular flux.
  TEST(out.write_angular_flux(state, quad));
  cout << "...wrote angular flux." << endl;

  out.finalize();
  cout << "...finalized out." << endl;

  return 0;
}

// Test writing material into different directories
int test_SiloOutput_directory(int argc, char *argv[])
{
  // Create mesh
  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0);
  cm[1] = 10;
  vec_int mt(1, 0);
  Mesh2D::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));

  // Create field of fata
  vec_dbl data(mesh->number_cells(), 1.234);

  // Create the SiloOutput.
  SiloOutput out(mesh);

  // Open the file, and write the mesh.
  TEST(out.initialize("test_dir.silo"));

  // Write out data
  TEST(out.write_scalar_field("SCALARFIELD", data));

  // Create a new directory and go to it
  TEST(out.make_directory("subdirectory"));
  TEST(out.set_directory("subdirectory"));

  // Write out the data again
  vec_dbl data2(data);
  TEST(out.write_scalar_field("SUBDIRFIELD2", data2))

  out.finalize();

  return 0;
}
//---------------------------------------------------------------------------//
//              end of test_SiloOutput.cc
//---------------------------------------------------------------------------//
