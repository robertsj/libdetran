//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utilities/DBC.hh
 * \brief  Design-by-Contract macros and abstract Object class.
 * \author Jeremy Roberts
 */

#ifndef DBC_HH
#define DBC_HH

#include "GenException.hh"
#include "detran_config.hh"

#include <iostream>
#include <string>

namespace detran_utilities
{

// Object class is the base class for all
// objects in the system. All classes inheriting from this class need 
// to define a method IsValid. This method should perform a
// consistency check on the state of the object. Note that 
// this method needs to be defined only when a debug build is made

//===========================================================================//
/*!
 * \class Object
 * \brief Abstract class for all objects following DBC
 *
 * This class and the associated macros are largely based on the nice 
 * tutorial given at:
 * http://eventhelix.com/realtimemantra/object_oriented/design_by_contract.htm 
 *
 */
//===========================================================================//
class Object
{

public:

  /// Virtual destructor
  ~Object(){}

#ifdef DETRAN_ENABLE_DEBUG
  /// Tests the object's state.
  virtual bool is_valid() const = 0;
#endif

};

//===========================================================================//
/*!
 * \page  DBC
 * \brief Design-By-Contract
 *
 * DBC is a framework for ensuring consistency between a client's request
 * and a server's response.  Here, we use a few macro's that help
 * implement the idea:
 *  - Require: precondition upon entering a method that must be satisfied;
 *             if not satisfied, the client is to blame
 *  - Ensure:  postcondition that must be satisfied for the server to be
 *             correct
 *  - IsValid: check an object's state for correctness
 *  - Assert:  a general assertion macro
 *  - Insist:  Conditions that must always be true (i.e. these are always
 *             checked)
 *
 * See the nice tutorial at:
 * http://eventhelix.com/realtimemantra/object_oriented/design_by_contract.htm
 *
 */
//===========================================================================//

#ifdef DETRAN_ENABLE_DEBUG

// DBC Macros
#define Assert(c)     if (!(c)) throw detran_utilities::GenException( __LINE__, __FILE__,#c)
#define Require(c)    Assert(c)
#define Ensure(c)     Assert(c)

// DBC Macros with additional message
#define Assert_msg(c, m)     if (!(c)) throw detran_utilities::GenException( __LINE__, __FILE__,std::string(#c)+" "+std::string(m))
#define Require_msg(c, m)    Assert_msg(c, m)
#define Ensure_msg(c, m)     Assert_msg(c, m)

#else

#define Assert(c)   ((void) 0)
#define Require(c)  ((void) 0)
#define Ensure(c)   ((void) 0)
#define Assert_msg(c, m)   ((void) 0)
#define Require_msg(c, m)  ((void) 0)
#define Ensure_msg(c, m)   ((void) 0)

#endif

#define Insist(c,m)   if (!(c)) {std::cerr << m << std::endl; throw detran_utilities::GenException( __LINE__, __FILE__,#c);}

} // end namespace detran_utilities

#endif // DBC_HH


//---------------------------------------------------------------------------//
//              end of DBC.hh
//---------------------------------------------------------------------------//
