//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Iterators.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_Iterators_Reversible)

#include "TestDriver.hh"
#include "utilities/Iterators.hh"
#include <vector>
#include <cstdio>

using namespace detran_test;
using namespace detran_utilities;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_Iterators_Reversible(int argc, char *argv[])//not testing anything
{

  std::vector<int>  v = {1, 2, 3};

  Reversible<std::vector<int>> forward(&v, true);

  int i;
  for (i = 0; i < 3; i++, forward++)
  {
    std::printf(" vf = %i\n", *forward);
    TEST(*forward == v[i]);
  }
  i = 0;
  for (auto it = forward.begin(); it != forward.end(); ++it, ++i)
  {
    std::printf(" vf = %i   %i\n", *it, v[i]);
    TEST(*it == v[i]);
  }
  i = 0;
  for (auto &it: forward)
  {
    TEST(it == v[i]);
    ++i;
  }

  Reversible<std::vector<int>> backward(&v, false);
  for (i = 0; i < 3; i++, backward++)
  {
    std::printf(" vf = %i\n", *backward);
    TEST(*backward == v[2 - i]);
  }
  i = 0;
  for (auto it = backward.begin(); it != backward.end(); ++it, ++i)
  {
    std::printf(" vf = %i   %i\n", *it, v[2 - i]);
    TEST(*it == v[2 - i]);
  }
  i = 0;
  for (auto &it : backward)
  {
    TEST(it == v[2 - i]);
    ++i;
  }


  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Iterators.cc
//----------------------------------------------------------------------------//
