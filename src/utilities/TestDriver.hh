//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TestDriver.hh
 * \author Jeremy Roberts
 * \date   Jul 12, 2011
 * \brief  Simple functions and macro for testing and printing failures.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef TESTDRIVER_HH_
#define TESTDRIVER_HH_

// Utilities
#include "DBC.hh"
#include "SoftEquivalence.hh"

// System
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>


/*!
 *  \page testing Detran Testing Typedefs and Macros
 *
 *  \section testcreation Creating a Unit Test Set
 *
 *  To produce a set of unit tests,
 *  the client must define in the main test file
 *  a list of the form
 *  \code
    #define TEST_LIST  \
        FUNC(test_one) \
        FUNC(test_two)
    #include "TestDriver.hh"
    #include "MyClassBeingTested.hh"
    #include "AnyNeededHelperClasses.hh"
    using namespace detran_test;
    using namespace detran;
    using namespace std;
    int main(int argc, char *argv[])
    {
      RUN(argc, argv);
    }
    int test_one() // and so on
    \endcode
 *
 *  Note, the TestDriver header must be included *after* the
 *  list is given.  There
 *  are ways around this, but if I'm to spend any more time on this, I'd
 *  be compelled to switch to someone else's unit test framework.
 *
 *  Once a set of tests is defined in this way, CTest can be used.  For
 *  example, in the CMakeLists.txt file associated with the class
 *  being tested, we would have entries like
 *
 *  \code
    ADD_EXECUTABLE(test_MyClass test_MyClass.cc)
    TARGET_LINK_LIBRARIES(test_MyClass some_library_myclass_needs)
    ADD_TEST( test_MyClass_one  test_MyClass 0)
    ADD_TEST( test_MyClass_two  test_MyClass 1)
 *  \endcode
 *
 *  where each test function in the source file is given its own
 *  CTest name and an integer index.  From there, various
 *  approaches to using CTest are possible, from one-off tests
 *  to dedicated nightly regression tests.  A full test set
 *  is under development.
 *
 *  \section defineunittest Defining a Unit Test
 *
 *  Above, we showed how to produce the driver file for one or
 *  more unit tests.  The actual unit test takes the form
 *
 *  \code
    int test_one()
    {
      // Do some setup.
      TEST(some_boolean_expression);
      // Do something else.
      TEST(some_other_boolean_expression);
      // Make sure to return 0 if all is well.
      return 0;
    }
 *  \endcode
 *
 *  Here, the TEST macro returns from the test function
 *  with a nonzero value, indicating to the driver
 *  function run that a specific test function has
 *  failed.  The rest of the functions, if any, will
 *  still be completed.
 */

/// \{

/// Typedefs for test arrays.
typedef int (*test_pointer)(void);
typedef std::string test_name;

/// Create test list.
#define FUNC( _name ) int _name(void);
TEST_LIST
#undef FUNC

/// Convert test list to array of function pointers.
#define FUNC( _name ) &_name,
test_pointer test_table[] =
{
  TEST_LIST
};
#undef FUNC

/// Convert test list to array of names of functions.
#define FUNC( _name ) #_name,
test_name test_names[] =
{
  TEST_LIST
};
#undef FUNC

/// Compute the number of tests in the list.
#define NUMBER_TESTS() ((sizeof(test_table)/sizeof(0[test_table])) \
                       / ((size_t)(!(sizeof(test_table)            \
                       % sizeof(0[test_table])))))

/// \}

namespace detran_test
{

/// \name Utility functions
/// \{

/// integer to string
static std::string itoa(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

/// double to string
static std::string dtoa(double d)
{
  std::stringstream out;
  out << d;
  return out.str();
}

/// \}

class TestDriver
{

public:

  /*!
   *  \brief Run the test specified at the command line.
   *
   *  The user test_XYZ.cc defines the test function list, and
   *  the command line contains the test indentifier.
   */
  static int run(int argc, char *argv[])
  {

    // Get test number
    Insist(argc==2, "Tests should have only one argument!");
    d_test = std::atoi(argv[1]);

    // Check if valid test identifier.
    Insist(d_test >= 0, "Test id must be nonnegative!");
    Insist(d_test < NUMBER_TESTS(),
           "Test id must be less than the number of tests defined!");

    // Announce test and return the identifier.
    std::cout << "running: " << test_names[d_test]
              << " ......" << std::endl;

    // Evaluate the test, and catch any exceptions; we don't want
    // all the tests to crash (if there are more than one in some potential
    // use case).
    try
    {
      std::cout << "...... " << test_names[d_test];
      if (!(*test_table[d_test])())
      {
        std::cout << " passed ";
      }
      else
      {
        std::cout << " failed ";
      }
      std::cout << std::endl;
    }
    catch (std::exception &err)
    {
      std::cout << std::endl;
      std::cout << "ERROR: While testing " << test_names[d_test]
           << " " << err.what()
           << std::endl;
      return 1;
    }
    catch ( ... )
    {
      return 1;
      std::cout << "ERROR: While testing " << test_names[d_test]
           << "An UNKNOWN exception was thrown."
           << std::endl;
      return 2;
    }
    return 1;
  }

  static int number_fails()
  {
    return d_number_fails;
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
        << "*** failed!         " << std::endl << std::endl;
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
    {
      return false;
    }
    else
    {
      return fail(cond, file, func, line);
    }
  }

  static int d_number_tests;
  static int d_number_fails;
  static int d_test;

};

int TestDriver::d_number_tests = 0;
int TestDriver::d_number_fails = 0;
int TestDriver::d_test = 0;


} // end namespace testing

/// Evaluate a statement expected to be true.
#define TEST(c)      if(TestDriver::                                 \
 test( c, #c, __FILE__, __func__, __LINE__ )) return 1;

/// Evaluate a statement expected to be false.
#define TESTFALSE(c) if(!TestDriver::                                 \
 test( c, #c, __FILE__, __func__, __LINE__ )) return 1;

/// The single function call needed in test_XYX.cc.
#define RUN(argv, argc) return TestDriver::run(argv, argc);

#endif /* TESTING_HH_ */

//---------------------------------------------------------------------------//
//              end of Testing.hh
//---------------------------------------------------------------------------//
