//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_InputDB.cc
 * \author Jeremy Roberts
 * \date   Jul 14, 2011
 * \brief  Test of InputDB
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST           \
        FUNC(test_InputDB)

// Detran headers
#include "TestDriver.hh"
#include "InputDB.hh"

// System headers
#include <iostream>
#include <cstdlib>

using namespace detran_utils;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

// Test definitions.

int test_InputDB()
{
 InputDB::SP_input db;
 db = new InputDB();
 db->put<int>("number_groups", 2);
 TEST(db->check("number_groups"));
 int ng = db->get<int>("number_groups");
 TEST(ng == 2);
 TEST(db->check("i_am_not_here"));
 return 0;
}


//---------------------------------------------------------------------------//
//              end of testTesting.cc
//---------------------------------------------------------------------------//
