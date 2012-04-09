//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   input_fixture.hh
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  An input for testing.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef INPUT_FIXTURE_HH_
#define INPUT_FIXTURE_HH_

// Detran includes
#include "InputDB.hh"

// System includes

namespace detran_test
{

typedef detran::InputDB::SP_input SP_input;

/*!
 *  \brief Create a common input for transport test cases.
 *
 *  This contains several default values that can be changed
 *  where needed.
 */
static SP_input input_fixture()
{
  // Create the new database.
  SP_input input;
  input = new detran::InputDB();

  // Enter some basic things.
  input->put<int>("number_groups", 2);
  input->put<int>("dimension", 1);

  // Ensure a valid input.
  Ensure(input->is_valid());

  // Return the fixture.
  return input;
}

} // end namespace detran_test

#endif /* INPUT_FIXTURE_HH_ */

//---------------------------------------------------------------------------//
//              end of input_fixture.hh
//---------------------------------------------------------------------------//
