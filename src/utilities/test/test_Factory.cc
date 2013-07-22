//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Factory.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST          \
        FUNC(test_Factory)

#include "TestDriver.hh"
#include "utilities/Factory.hh"
#include <vector>

using namespace std;
using namespace detran_test;
using namespace detran_utilities;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

/// The base class
class Foo
{
public:
  typedef SP<Foo> SP_foo;
  // A REQUIRED typedef that tells the factory what the function type is
  typedef SP_foo (*CreateFunction)(int n);
  // A REQUIRED typedef that defines the factory
  typedef Factory<Foo>  Factory_T;
  Foo(int n) : d_n(n) { }
  virtual ~Foo(){}
  virtual int value() = 0;
  /**
   *  A SUGGESTED method that serves as the client's point of access. This
   *  can be in the base class or outside, as illustrated below
   */
  static SP_foo Create(const std::string &key, int n)
  {
    return (Factory_T::Instance().GetCreateFunction(key))(n);
  }
protected:
  int d_n;
};

// Example of non-member SUGGESTED method for the client's point of access
inline Foo::SP_foo Create(const std::string &key, int n)
{
  return (Foo::Factory_T::Instance().GetCreateFunction(key))(n);
}

/// The creation function template with user-defini
template <typename D>
Foo::SP_foo Create(int n)
{
  return Foo::SP_foo(new D(n));
}


class Moo: public Foo
{
public:
  Moo(int n) : Foo(n) {}
  virtual int value() { return d_n;}
};
REGISTER_CLASS(Foo, Moo, "Moo")

class Noo: public Foo
{
public:
  Noo(int n) : Foo(n) {}
  virtual int value() { return 2 * d_n;}
};
REGISTER_CLASS(Foo, Noo, "Noo")

#define PRINT(c) cout << c << endl;
int test_Factory(int argc, char *argv[])
{
  // illustrates creation as a base method (probably a preferred approach)
  Foo::SP_foo moo = Foo::Create("Moo", 3);
  PRINT("hello")
  TEST(moo);
  // illustrates creation as a standalone factory method
  Foo::SP_foo noo = Create("Noo", 3);
  TEST(noo);

  std::cout << "moo=" << moo->value() << std::endl;
  std::cout << "noo=" << noo->value() << std::endl;

  TEST(moo->value() == 3)
  TEST(noo->value() == 6)

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Factory.cc
//----------------------------------------------------------------------------//
