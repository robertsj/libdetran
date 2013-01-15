//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_Vector.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of Vector class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_Vector)       \
        FUNC(test_Vector_resize)

#include "TestDriver.hh"
#include "callow/vector/Vector.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/Definitions.hh"
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
  {
    Vector v(10, 1.0);
    for (int i = 0; i < 10; i++)
    {
      TEST(soft_equiv(v[i], 1.0));
    }
    v[5] = 5.0;
    TEST(soft_equiv(v[5], 5.0));
    v[5] = 1.0;

    Vector y(10, 2.0);
    double val = v.dot(y);
    double val2 = 0.0;
    for (int i = 0; i < 10; i++)
    {
      val2 += v[i]*y[i];
    }
    cout << val2 << endl;
    TEST(soft_equiv(val, 20.0));

    v.set(1.23);
    TEST(soft_equiv(v[0], 1.23));
    v.scale(2.0);
    TEST(soft_equiv(v[0], 2.46));
  }

  {
    detran_utilities::vec_dbl v(5, 1.0);
    Vector V(v);
    TEST(V.size() == 5);
    TEST(soft_equiv(V[4], 1.0));
  }

  {
    double v[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    Vector V(5, v);
    TEST(V.size() == 5);
    TEST(soft_equiv(V[4], 1.0));
    double v_L1 = V.norm(L1);
    TEST(soft_equiv(v_L1, 5.0));
    double v_L2 = V.norm(L2);
    TEST(soft_equiv(v_L2, 2.236067977499790));
    double v_LI = V.norm(LINF);
    TEST(soft_equiv(v_LI, 1.0));
  }

  return 0;
}

int test_Vector_resize(int argc, char *argv[])
{
  typedef Vector Vec_T;
  Vec_T v;
  v.resize(5, 1.0);
  v.set(3.0);
  v.display();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Vector.cc
//---------------------------------------------------------------------------//
