//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  test_InputDB.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "InputDB.hh"
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace detran_utilities;
using namespace std;

TEST(InputDB, Basic)
{
 InputDB::SP_input db(new InputDB());

 // Check that when we insert, it shows up, and we can get it.

 // int
 db->put<int>("number_groups", 2);
 EXPECT_TRUE(db->check("number_groups"));
 int ng = db->get<int>("number_groups");
 EXPECT_EQ(ng, 2);

 // double
 db->put<double>("inner_tolerance", 0.0001);
 EXPECT_TRUE(db->check("inner_tolerance"));
 double tol = db->get<double>("inner_tolerance");
 EXPECT_NEAR(tol, 0.0001, 1.0e-12);

 // vecint
 vec_int vint(10, 1);
 db->put<vec_int>("vector_of_ints", vint);
 EXPECT_TRUE(db->check("vector_of_ints"));
 vec_int vint2 = db->get<vec_int>("vector_of_ints");
 EXPECT_EQ(vint2.size(), 10);
 for (int i = 0; i < 10; i++)
 {
    EXPECT_EQ(vint[i], vint2[i]);
 }

 // vecdbl
 vec_dbl vdbl(10, 1.0);
 db->put<vec_dbl>("vector_of_dbls", vdbl);
 EXPECT_TRUE(db->check("vector_of_dbls"));
 vec_dbl vdbl2 = db->get<vec_dbl>("vector_of_dbls");
 EXPECT_EQ(vdbl2.size(), 10);
 for (int i = 0; i < 10; i++)
 {
   EXPECT_NEAR(vdbl[i], vdbl2[i] ,1.0e-12);
 }

 // string
 db->put<string>("equation", "dd");
 EXPECT_TRUE(db->check("equation"));
 string eq = db->get<string>("equation");
 EXPECT_EQ(eq, "dd");

 // input db
 InputDB::SP_input db2(new InputDB("DB2"));
 db2->put<int>("testint", 99);
 db2->put<double>("testdbl", 1.23);
 db->put<InputDB::SP_input>("db2", db2);
 EXPECT_TRUE(db->check("db2"));
 InputDB::SP_input db2check;
 db2check = db->get<InputDB::SP_input>("db2");
 EXPECT_TRUE(db2check != NULL);
 int testint = db2check->get<int>("testint");
 EXPECT_EQ(testint, 99);
 double testdbl = db2check->get<double>("testdbl");
 EXPECT_NEAR(testdbl, 1.23, 1.0e-12);

 // Test that something is not there.
 EXPECT_FALSE(db->check("i_am_not_here"));
}

TEST(InputDB, Serialize)
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
    EXPECT_TRUE(db->check("number_groups"));
    EXPECT_EQ(db->get<int>("number_groups"), 2);
    EXPECT_TRUE(db->check("inner_tolerance"));
    EXPECT_NEAR(db->get<double>("inner_tolerance"), 0.0001, 1.0e-12);
    EXPECT_TRUE(db->check("vector_of_ints"));
    vec_int vi = db->get<vec_int>("vector_of_ints");
    EXPECT_EQ(vi.size() == 10);
    for (int i = 0; i < 10; i++)
      EXPECT_EQ(vi[i], 1);
    EXPECT_TRUE((db->check("vector_of_dbls"));
    vec_dbl vd = db->get<vec_dbl>("vector_of_dbls");
    EXPECT_EQ(vd.size(), 10);
    for (int i = 0; i < 10; i++)
      EXPECT_NEAR(vd[i], 1.0, 1.0e-12);
    EXPECT_TRUE(db->check("equation"));
    string eq = db->get<string>("equation");
    EXPECT_EQ(eq, db->get<string>("equation"));
    EXPECT_FALSE(db->check("i_am_not_here"));
  }
#endif
}

//---------------------------------------------------------------------------//
//              end of test_InputDB.cc
//---------------------------------------------------------------------------//
