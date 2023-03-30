//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test.cc
 *  @brief Test of Matrix class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "matrix_fixture.hh"
#include "matrix/Matrix.hh"
#include "matrix/SparseRow.hh"
#include "utils/Initialization.hh"
#include <iostream>
#include <algorithm>

using namespace callow;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
TEST(SparceRow, Basic)
{
  const int n = 5;
  Matrix::SP_matrix A = test_matrix_1(n);

  double l = -0.20;
  double d =  0.50;
  double u = -0.21;

  {
    SparseRow::element_t element(1, 2.0);
    EXPECT_EQ(element.first, 1);
    EXPECT_EQ(element.second, 2.0);
  }

  {
	//              0  1  2  3  4  5  6  7  8
    double a[10] = {1, 0, 3, 0, 2, 0, 6, 4, 5, };
    SparseRow row(a, 10);
    SparseRow row2 = row;
    std::sort(row.begin(), row.end(), SparseRow::compare_v);
    EXPECT_EQ(row.begin()->first, 6);
    EXPECT_EQ((row.begin()+5)->first, 0);
    EXPECT_EQ((row.begin()+5)->second, 1);
    std::sort(row2.range(0, 5).first, row2.range(0, 5).second, SparseRow::compare_v);
    EXPECT_EQ(row2.begin()->first, 2);
    EXPECT_EQ((row2.begin()+5)->first, 8);
    EXPECT_EQ((row2.begin()+5)->second, 5);
  }

  {
    // test row from sparse matrix and copy
    A->display();
    SparseRow row(*A, 1);
    SparseRow wow = row;
    row.display();
    EXPECT_EQ(row.size(), n);
    EXPECT_EQ(row.number_nonzeros(), 3);
    EXPECT_NEAR(row.elements()[0].second, l, 1.0e-12);
    EXPECT_NEAR(row.elements()[1].second, d, 1.0e-12);
    EXPECT_NEAR(row.elements()[2].second, u, 1.0e-12);
    EXPECT_EQ(row.elements()[0].first, 0);
    EXPECT_EQ(row.elements()[1].first, 1);
    EXPECT_EQ(row.elements()[2].first, 2);
    row.clear();
    EXPECT_EQ(row.number_nonzeros(), 0);
    EXPECT_EQ(wow.size(), n);
    EXPECT_EQ(wow.number_nonzeros(), 3);
    EXPECT_NEAR(wow.elements()[0].second, l, 1.0e-12);
    EXPECT_NEAR(wow.elements()[1].second, d, 1.0e-12);
    EXPECT_NEAR(wow.elements()[2].second, u, 1.0e-12);
    EXPECT_EQ(wow.elements()[0].first, 0);
    EXPECT_EQ(wow.elements()[1].first, 1);
    EXPECT_EQ(wow.elements()[2].first, 2);
    // test single insert
    wow.display();
    wow.insert(0, 1.0, SparseRow::ADD);
    wow.display();
    EXPECT_NEAR(wow.elements()[0].second, 1.0+l, 1.0e-12);
    EXPECT_TRUE(wow.insert(4, 99.0, SparseRow::INSERT));
    EXPECT_TRUE(wow.verify());
    wow.display();
    EXPECT_NEAR(wow.elements()[3].second, 99.0, 1.0e-12);
    EXPECT_EQ(wow.number_nonzeros(), 4);
    wow.clear();
    wow.insert(3, 1.0, SparseRow::INSERT);
  }

  {
	// test row from array
    double arr0[n] = {1, 0, 2, 0, 0};
    double arr1[n] = {0, 3, 4, 0, 5};
    SparseRow row0(arr0, n);
    SparseRow row1(arr1, n);
    EXPECT_EQ(row0.number_nonzeros(), 2);
    EXPECT_EQ(row1.number_nonzeros(), 3);
    EXPECT_NEAR(row0.elements()[0].second, 1., 1.0e-12);
    EXPECT_NEAR(row0.elements()[1].second, 2., 1.0e-12);
    EXPECT_NEAR(row1.elements()[0].second, 3., 1.0e-12);
    EXPECT_NEAR(row1.elements()[1].second, 4., 1.0e-12);
    EXPECT_NEAR(row1.elements()[2].second, 5., 1.0e-12);
    EXPECT_EQ(row0.elements()[0].first, 0);
    EXPECT_EQ(row0.elements()[1].first, 2);
    EXPECT_EQ(row1.elements()[0].first, 1);
    EXPECT_EQ(row1.elements()[1].first, 2);
    EXPECT_EQ(row1.elements()[2].first, 4);

    // test full addition [x,x,x,0,x]
    SparseRow row2(row0);                         // [1, 0, 2, 0, 0]
    row2.insert(0, n, row1, 1, SparseRow::ADD);   // [1, 3, 6, 0, 5]
    row2.display();
    EXPECT_EQ(row2.number_nonzeros(), 4);
    EXPECT_NEAR(row2.elements()[0].second, 1., 1.0e-12);
    EXPECT_NEAR(row2.elements()[1].second, 3., 1.0e-12);
    EXPECT_NEAR(row2.elements()[2].second, 6., 1.0e-12);
    EXPECT_NEAR(row2.elements()[3].second, 5., 1.0e-12);
    EXPECT_EQ(row2.elements()[0].first, 0);
    EXPECT_EQ(row2.elements()[1].first, 1);
    EXPECT_EQ(row2.elements()[2].first, 2);
    EXPECT_EQ(row2.elements()[3].first, 4);


    // test partial addition to do row2[a:b]
    row2.insert(0, 2, row0, -1, SparseRow::ADD);  // [0, 3, 6, 0, 5]
    for (auto it = row0.range(0, 2).first; it < row0.range(0, 2).second; ++it)
      cout << it->first << " " << it->second << endl;


    row2.display();
    row2.update();
    EXPECT_EQ(row2.number_nonzeros(), 3);
    EXPECT_NEAR(row2.elements()[0].second, 3., 1.0e-12);
    EXPECT_NEAR(row2.elements()[1].second, 6., 1.0e-12);
    EXPECT_NEAR(row2.elements()[2].second, 5., 1.0e-12);
    EXPECT_EQ(row2.elements()[0].first, 1);
    EXPECT_EQ(row2.elements()[1].first, 2);
    EXPECT_EQ(row2.elements()[2].first, 4);
  }

  {
    // test iterator
    double arr1[n] = {0, 3, 4, 0, 5};
    SparseRow row1(arr1, n);

    for (SparseRow::iterator it = row1.begin();  it != row1.end(); ++it)
    {
      cout << it->first << " " << it->second << endl;
    }
  }
  
  callow_finalize();
}


//----------------------------------------------------------------------------//
//              end of test_SparseRow.cc
//----------------------------------------------------------------------------//
