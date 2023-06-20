//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Iterators.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_Iterators_Reversible)

#include "utilities/Iterators.hh"
#include <vector>
#include <cstdio>

using namespace detran_utilities;


//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

TEST(Iterators, Reversible)
{
  std::vector<int>  v = {1, 2, 3};

  Reversible<std::vector<int>> forward(&v, true);

  int i;
  for (i = 0; i < 3; i++, forward++)
  {
    std::printf(" vf = %i\n", *forward);
    EXPECT_EQ(*forward, v[i]);
  }
  i = 0;
  for (auto it = forward.begin(); it != forward.end(); ++it, ++i)
  {
    std::printf(" vf = %i   %i\n", *it, v[i]);
    EXPECT_EQ(*it, v[i]);
  }
  i = 0;
  for (auto &it: forward)
  {
    EXPECT_EQ(it, v[i]);
    ++i;
  }

  Reversible<std::vector<int>> backward(&v, false);
  for (i = 0; i < 3; i++, backward++)
  {
    std::printf(" vf = %i\n", *backward);
    EXPECT_EQ(*backward, v[2 - i]);
  }
  i = 0;
  for (auto it = backward.begin(); it != backward.end(); ++it, ++i)
  {
    std::printf(" vf = %i   %i\n", *it, v[2 - i]);
    EXPECT_EQ(*it, v[2 - i]);
  }
  i = 0;
  for (auto &it : backward)
  {
    EXPECT_EQ(it, v[2 - i]);
    ++i;
  }
}

//----------------------------------------------------------------------------//
//              end of test_Iterators.cc
//----------------------------------------------------------------------------//
