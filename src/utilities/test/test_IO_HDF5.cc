//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_IO_HDF5.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of IO_HDF5 class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_IO_HDF5)

// Detran headers
#include "TestDriver.hh"
#include "IO_HDF5.hh"

// Setup
/* ... */

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_IO_HDF5(int argc, char *argv[])
{

  // Create the input.
  InputDB::SP_input input(new InputDB());

  // Put some integers
  input->put<int>("int0", 0);
  input->put<int>("int1", 1);

  // Put some doubles
  input->put<double>("double0_123", 0.123);
  input->put<double>("double0_234", 0.234);

  // Put some strings
  input->put<string>("string_test", "test");
  input->put<string>("hdf5_input",  "lala.h5");

  // Put vectors.
  vec_int v_int(4, 0);
  v_int[1] = 1;
  v_int[2] = 2;
  v_int[3] = 3;
  vec_dbl v_dbl(4, 0.0);
  v_dbl[1] = 2.0;
  v_dbl[2] = 4.0;
  v_dbl[3] = 6.0;
  input->put<vec_int>("vec_int_0123", v_int);
  input->put<vec_dbl>("vec_dbl_0246", v_dbl);

  // Create an IO_HDF5
  IO_HDF5 io("test.h5");

  // Write to the HDF5 file.  The input has the filename
  // to write out, or a default is used.
  io.write(input);
  io.close();

  // Create a new input database and read in from the
  InputDB::SP_input input2(new InputDB());
  io.read(input2, "lala.h5");

  // Verify that the input = input2
  TEST(input2->get<int>("int0") == 0);
  TEST(input2->get<int>("int1") == 1);


  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_IO_HDF5.cc
//---------------------------------------------------------------------------//
