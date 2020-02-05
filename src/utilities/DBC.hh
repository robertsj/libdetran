//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   utilities/DBC.hh
 *  @brief  Design-by-Contract macros and abstract Object class.
 *  @author Jeremy Roberts
 */

#ifndef detran_utilities_DBC_HH
#define detran_utilities_DBC_HH

#include "GenException.hh"
#include "detran_config.hh"

#include <iostream>
#include <string>
#include <typeinfo>
#include <sstream>

namespace detran_utilities
{

//===========================================================================//
/**
 *  @page  DBC
 *  @brief Design-By-Contract
 *
 *  DBC is a framework for ensuring consistency between a client's request
 *  and a server's response.  Here, we use a few macro's that help
 *  implement the idea:
 *   - Require: precondition upon entering a method that must be satisfied;
 *              if not satisfied, the client is to blame
 *   - Ensure:  postcondition that must be satisfied for the server to be
 *              correct
 *   - Assert:  a general assertion macro
 *   - Insist:  Conditions that must always be true (i.e. these are always
 *              checked)
 *
 *  See the nice tutorial at:
 *  http://eventhelix.com/realtimemantra/object_oriented/design_by_contract.htm
 *
 */
//===========================================================================//

/// Compile time boolean (straight from Alexandrescu's Modern C++ Design
template <bool> struct StaticChecker {StaticChecker(...);};
template <> struct StaticChecker<false> { };

#ifdef DETRAN_ENABLE_DEBUG

// DBC Macros
#define Assert(c) if (!(c)) throw detran_utilities::GenException( __LINE__, __FILE__,#c)
#define Require(c) Assert(c)
#define Ensure(c)  Assert(c)

// Verbose DBC Macros (i.e. with an additional message)
#define Assertv(c, m) if (!(c)) throw detran_utilities::GenException( __LINE__, __FILE__, std::string(#c) +", " + std::string(m))
#define Requirev(c, m) Assertv(c, m)
#define Ensurev(c, m)  Assertv(c, m)

// Macro for compile time assertions (from Alexandrescu)
#define StaticAssert(c)                                                        \
{                                                                              \
  class ERROR_Compile_time{};                                                  \
  (void)sizeof detran_utilities::StaticChecker<(c)!=0>((ERROR_Compile_time()));\
}
#define StaticAssertv(c, m)                                             \
{                                                                       \
  class ERROR_##m{};                                                    \
  (void)sizeof detran_utilities::StaticChecker<(c)!=0>((ERROR_##m()) ); \
}

#else

#define Assert(c)           ((void) 0);
#define Require(c)          ((void) 0);
#define Ensure(c)           ((void) 0);
#define StaticAssert(c)     ((void) 0);
#define Assertv(c, m)       ((void) 0);
#define Requirev(c, m)      ((void) 0);
#define Ensurev(c, m)       ((void) 0);
#define StaticAssertv(c, m) ((void) 0);

#endif

#define Insist(c,m) if (!(c)) {std::cerr << m << std::endl; throw detran_utilities::GenException( __LINE__, __FILE__, #c);}

template <class T>
inline std::string as_string(T v)
{
  std::ostringstream o;
  if (!(o << v))
    THROW("Error converting to string from " + std::string(typeid(T).name()))
  return o.str();
}

#define AsString(c) detran_utilities::as_string(c)

} // end namespace detran_utilities

#endif // detran_utilities_DBC_HH

//---------------------------------------------------------------------------//
//              end of DBC.hh
//---------------------------------------------------------------------------//
