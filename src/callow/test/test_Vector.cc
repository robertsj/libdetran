//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Vector.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of Vector class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_Vector)

#include "TestDriver.hh"
#include "vector/Vector.hh"
#include "utils/Initialization.hh"
#include <iostream>

using namespace callow;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// Test of basic public interface
int test_Vector(int argc, char *argv[])
{
  typedef Vector<double> vec_dbl;

  vec_dbl v(10, 1.0);
  for (int i = 0; i < 10; i++)
  {
    TEST(soft_equiv(v[i], 1.0));
  }
  v[5] = 5.0;
  TEST(soft_equiv(v[5], 5.0));
  v[5] = 1.0;

  vec_dbl y(10, 2.0);
  double val = v.dot(y);
  double val2 = 0.0;
  for (int i = 0; i < 10; i++)
  {
    val2 += v[i]*y[i];
  }
  cout << val2 << endl;
  TEST(soft_equiv(val, 20.0));



  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Vector.cc
//---------------------------------------------------------------------------//
