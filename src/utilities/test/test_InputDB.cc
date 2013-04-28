//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  test_InputDB.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_InputDB)           \
        FUNC(test_InputDB_serialize)

// Detran headers
#include "TestDriver.hh"
#include "InputDB.hh"

// System headers
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace detran_test;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

// Test definitions.

int test_InputDB(int argc, char *argv[])
{
 InputDB::SP_input db;
 db = new InputDB();

 // Check that when we insert, it shows up, and we can get it.

 // int
 db->put<int>("number_groups", 2);
 TEST(db->check("number_groups"));
 int ng = db->get<int>("number_groups");
 TEST(ng == 2);

 // double
 db->put<double>("inner_tolerance", 0.0001);
 TEST(db->check("inner_tolerance"));
 double tol = db->get<double>("inner_tolerance");
 TEST(soft_equiv(tol, 0.0001));

 // vecint
 vec_int vint(10, 1);
 db->put<vec_int>("vector_of_ints", vint);
 TEST(db->check("vector_of_ints"));
 vec_int vint2 = db->get<vec_int>("vector_of_ints");
 TEST(vint2.size() == 10);
 for (int i = 0; i < 10; i++)
 {
   TEST(vint[i] == vint2[i]);
 }

 // vecdbl
 vec_dbl vdbl(10, 1.0);
 db->put<vec_dbl>("vector_of_dbls", vdbl);
 TEST(db->check("vector_of_dbls"));
 vec_dbl vdbl2 = db->get<vec_dbl>("vector_of_dbls");
 TEST(vdbl2.size() == 10);
 for (int i = 0; i < 10; i++)
 {
   TEST(soft_equiv(vdbl[i], vdbl2[i]));
 }

 // string
 db->put<string>("equation", "dd");
 TEST(db->check("equation"));
 string eq = db->get<string>("equation");
 TEST(eq=="dd");

 // input db
 InputDB::SP_input db2(new InputDB("DB2"));
 db2->put<int>("testint", 99);
 db2->put<double>("testdbl", 1.23);
 db->put<InputDB::SP_input>("db2", db2);
 TEST(db->check("db2"));
 InputDB::SP_input db2check;
 db2check = db->get<InputDB::SP_input>("db2");
 TEST(db2check);
 int testint = db2check->get<int>("testint");
 TEST(testint == 99);
 double testdbl = db2check->get<double>("testdbl");
 TEST(soft_equiv(testdbl, 1.23));

 // Test that something is not there.
 TEST(!db->check("i_am_not_here"));

 return 0;
}

int test_InputDB_serialize(int argc, char *argv[])
{

#ifdef DETRAN_ENABLE_BOOST
  // Create and pack an input
  {

    // Create Input
    InputDB::SP_input db(new InputDB());
    db->put<int>("number_groups", 2);
    db->put<double>("inner_tolerance", 0.0001);
    vec_int vint(10, 1);
    db->put<vec_int>("vector_of_ints", vint);
    vec_dbl vdbl(10, 1.0);
    db->put<vec_dbl>("vector_of_dbls", vdbl);
    db->put<string>("equation", "dd");

    // Output
    std::ofstream ofs("test.inp");
    boost::archive::binary_oarchive oa(ofs);
    oa << db;
    ofs.close();
  }

  // Unpack and test an input
  {
    // Empty database
    InputDB::SP_input db;

    // Load
    std::ifstream ifs("test.inp");
    boost::archive::binary_iarchive ia(ifs);
    ia >> db;
    ifs.close();

    // Test
    TEST(db->check("number_groups"));
    TEST(db->get<int>("number_groups") == 2);
    TEST(db->check("inner_tolerance"));
    TEST(soft_equiv(db->get<double>("inner_tolerance"), 0.0001));
    TEST(db->check("vector_of_ints"));
    vec_int vi = db->get<vec_int>("vector_of_ints");
    TEST(vi.size() == 10);
    for (int i = 0; i < 10; i++)
      TEST(vi[i] == 1);
    TEST(db->check("vector_of_dbls"));
    vec_dbl vd = db->get<vec_dbl>("vector_of_dbls");
    TEST(vd.size() == 10);
    for (int i = 0; i < 10; i++)
      TEST(soft_equiv(vd[i], 1.0));
    TEST(db->check("equation"));
    string eq = db->get<string>("equation");
    TEST(eq == db->get<string>("equation"));
    TEST(!db->check("i_am_not_here"));
  }
#endif

 return 0;
}

//---------------------------------------------------------------------------//
//              end of test_InputDB.cc
//---------------------------------------------------------------------------//
