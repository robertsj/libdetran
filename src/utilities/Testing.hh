//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Testing.hh
 * \author Jeremy Roberts
 * \date   Jul 12, 2011
 * \brief  Simple functions and macro for testing and printing failures.
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef TESTING_HH_
#define TESTING_HH_

// System
#include <iostream>
#include <string>
#include <sstream>

// \todo Update test names to be more descriptive (e.g. what method they test)

namespace testing
{

/// integer to string
std::string itoa(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

/// double to string
std::string dtoa(double d)
{
  std::stringstream out;
  out << d;
  return out.str();
}


/*!
 *  \brief Prints a failure message.
 *
 *  \param cond  logical expression represented as a string
 *  \param file  file in which the test is called
 *  \param func  function in which the test is called
 *  \param line  source code line at which the test is called
 *  \return      true (indicating failure)
 *
 */
static bool fail(std::string const & cond,
                 std::string const & file,
                 std::string const & func,
                 int line)
{
  std::ostringstream message;
  message << "*** Test condition: " << cond << std::endl
          << "***        on line: " << itoa(line) << std::endl
          << "***    in function: " << func << std::endl
          << "***        in file: " << file << std::endl
          << "*** failed!         "         << std::endl
          << std::endl;
  std::cout << message.str();
  return true;
}

/*!
 *  \brief Tests an expression.
 *
 *  This function could be used by itself, but it is meant
 *  to be called via the macro \ref TEST, via
 *  \code
 *  // do stuff
 *  TEST(foo1==foo2); // returns 1 upon failure
 *  // do more stuff
 *  \endcode
 *  Using the macro helps eliminate some repetition in the
 *  tests, but a better canned approach (e.g. Google test)
 *  might be desirable.  In any case, these tests return
 *  to the main driver (w/ a printed message), and the
 *  main driver can track the various fails.  The main
 *  driver returns a nonzero if there is any failure. CTest
 *  works based on these nonzero values.
 *
 *  \param c     logical expression tested by user
 *  \param cond  same expression represented as a string
 *  \param file  file in which the test is called
 *  \param func  function in which the test is called
 *  \param line  source code line at which the test is called
 *  \return   true for pass or false for failure
 *
 */
static bool test(bool c,
                 std::string const & cond,
                 std::string const & file,
                 std::string const & func,
                 int line)
{
  if (c)
    return false;
  else
    return fail(cond, file, func, line);
}

} // end namespace testing

/// test macro
#define TEST(c) if(testing::test( c, #c, __FILE__, __func__, __LINE__ )) return 1;

/// print file name, which should describe test
#define TEST_ANNOUNCE() std::cout << "----" << __FILE__ << "----" << std::endl;

#endif /* TESTING_HH_ */

//---------------------------------------------------------------------------//
//              end of Testing.hh
//---------------------------------------------------------------------------//
