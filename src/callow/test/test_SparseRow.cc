//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test.cc
 *  @brief Test of Matrix class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST             \
        FUNC(test_SparseRow)

#include "TestDriver.hh"
#include "matrix_fixture.hh"
#include "matrix/Matrix.hh"
#include "matrix/SparseRow.hh"
#include "utils/Initialization.hh"
#include <iostream>
#include <algorithm>

using namespace callow;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;


int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
int test_SparseRow(int argc, char *argv[])
{
  callow_initialize(argc, argv);

  const int n = 5;
  Matrix::SP_matrix A = test_matrix_1(n);

  double l = -0.20;
  double d =  0.50;
  double u = -0.21;

  /**
  {
    // can we make a new vector from an iterator of another by reference?
    std::vector<int> v(5);
    for (int i = 0; i < 5; ++i)
      v[i] = i+1;
    auto it = v.begin()+2;
    std::vector<int> &y = *it;
  }
  */

  {
    SparseRow::element_t element(1, 2.0);
    TEST(element.first == 1);
    TEST(element.second == 2.0);
  }

  {
	//              0  1  2  3  4  5  6  7  8
    double a[10] = {1, 0, 3, 0, 2, 0, 6, 4, 5, };
    SparseRow row(a, 10);
    SparseRow row2 = row;
    row.display("before");
    std::sort(row.begin(), row.end(), SparseRow::compare);
    row.display("after");

    row2.display("before 2");
    std::sort(row2.range(0, 5).first, row2.range(0, 5).second, SparseRow::compare);
    row2.display("after 2");

  }

  callow_finalize();
  return 0;

  {
    // test row from sparse matrix and copy
    A->display();
    SparseRow row(*A, 1);
    SparseRow wow = row;
    row.display();
    TEST(row.size() == n);
    TEST(row.number_nonzeros() == 3);
    TEST(soft_equiv(row.elements()[0].second, l));
    TEST(soft_equiv(row.elements()[1].second, d));
    TEST(soft_equiv(row.elements()[2].second, u));
    TEST(row.elements()[0].first == 0);
    TEST(row.elements()[1].first == 1);
    TEST(row.elements()[2].first == 2);
    row.clear();
    TEST(row.number_nonzeros() == 0);
    TEST(wow.size() == n);
    TEST(wow.number_nonzeros() == 3);
    TEST(soft_equiv(wow.elements()[0].second, l));
    TEST(soft_equiv(wow.elements()[1].second, d));
    TEST(soft_equiv(wow.elements()[2].second, u));
    TEST(wow.elements()[0].first == 0);
    TEST(wow.elements()[1].first == 1);
    TEST(wow.elements()[2].first == 2);
    // test single insert
    wow.display();
    wow.insert(0, 1.0, SparseRow::ADD);
    wow.display();
    TEST(soft_equiv(wow.elements()[0].second, 1.0+l));
    TEST(wow.insert(4, 99.0, SparseRow::INSERT));
    TEST(wow.verify());
    wow.display();
    TEST(soft_equiv(wow.elements()[3].second, 99.0));
    TEST(wow.number_nonzeros() == 4);
    wow.clear();
    wow.insert(3, 1.0, SparseRow::INSERT);
  }

  // test row from array
  {
    double arr0[n] = {1, 0, 2, 0, 0};
    double arr1[n] = {0, 3, 4, 0, 5};
    SparseRow row0(arr0, n);
    SparseRow row1(arr1, n);
    TEST(row0.number_nonzeros() == 2);
    TEST(row1.number_nonzeros() == 3);
    TEST(soft_equiv(row0.elements()[0].second, 1.));
    TEST(soft_equiv(row0.elements()[1].second, 2.));
    TEST(soft_equiv(row1.elements()[0].second, 3.));
    TEST(soft_equiv(row1.elements()[1].second, 4.));
    TEST(soft_equiv(row1.elements()[2].second, 5.));
    TEST(row0.elements()[0].first == 0);
    TEST(row0.elements()[1].first == 2);
    TEST(row1.elements()[0].first == 1);
    TEST(row1.elements()[1].first == 2);
    TEST(row1.elements()[2].first == 4);

    // test full addition [x,x,x,0,x]
    SparseRow row2(row0);                         // [1, 0, 2, 0, 0]
    row2.insert(0, n, row1, 1, SparseRow::ADD);   // [1, 3, 6, 0, 5]
    row2.display();
    TEST(row2.number_nonzeros() == 4);
    TEST(soft_equiv(row2.elements()[0].second, 1.));
    TEST(soft_equiv(row2.elements()[1].second, 3.));
    TEST(soft_equiv(row2.elements()[2].second, 6.));
    TEST(soft_equiv(row2.elements()[3].second, 5.));
    TEST(row2.elements()[0].first == 0);
    TEST(row2.elements()[1].first == 1);
    TEST(row2.elements()[2].first == 2);
    TEST(row2.elements()[3].first == 4);


    // test partial addition to do row2[a:b]
    row2.insert(0, 2, row0, -1, SparseRow::ADD);  // [0, 3, 6, 0, 5]
    for (auto it = row0.range(0, 2).first; it < row0.range(0, 2).second; ++it)
      cout << it->first << " " << it->second << endl;


    row2.display();
    row2.update();
    TEST(row2.number_nonzeros() == 3);
    TEST(soft_equiv(row2.elements()[0].second, 3.));
    TEST(soft_equiv(row2.elements()[1].second, 6.));
    TEST(soft_equiv(row2.elements()[2].second, 5.));
    TEST(row2.elements()[0].first == 1);
    TEST(row2.elements()[1].first == 2);
    TEST(row2.elements()[2].first == 4);
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
  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_SparseRow.cc
//----------------------------------------------------------------------------//
