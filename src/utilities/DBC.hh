//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utilities/DBC.hh
 * \brief  Design-by-Contract macros and abstract Object class.
 * \author Jeremy Roberts
 */

#ifndef DBC_HH
#define DBC_HH

#include "detran_config.h"

namespace detran_utils
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

#ifdef DETRAN_ENABLE_DEBUG
    virtual bool is_valid() const = 0;
#endif

};



#ifdef DETRAN_ENABLE_DEBUG

#define ASSERT(c)   if (!(c)) abort_program(#c, __FILE__, __LINE__)
#define IS_VALID(c) ASSERT((obj) != NULL && (obj)->is_valid())
#define REQUIRE(c)  ASSERT(bool_expression)
#define ENSURE(c)   ASSERT(bool_expression)

#else

#define ASSERT(c)   ((void) 0)
#define IS_VALID(c) ((void) 0) 
#define REQUIRE(c)  ((void) 0)
#define ENSURE(c)   ((void) 0)

#endif

} // end namespace detran_utils

#endif // harness_DBC_hh


//---------------------------------------------------------------------------//
//              end of DBC.hh
//---------------------------------------------------------------------------//
