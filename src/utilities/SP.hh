//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utilities/SP.hh
 * \author Jeremy Roberts
 * \brief  Smart-Pointer (SP) class definition.
 * \note   Modified version of Tom Evan's SP class from Denovo.
 */
//---------------------------------------------------------------------------//

#ifndef SP_HH_
#define SP_HH_

#include "DBC.hh"

#ifdef DETRAN_ENABLE_BOOST
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#endif

#include <iostream>
#include <string>
#include <typeinfo>

namespace detran_utilities
{
 
//---------------------------------------------------------------------------//
/*!
 *  \class SPref
 *  \brief Reference counter for SP class.
 *
 *  This reference counter is thread safe, and allows SP's to be used
 *  as pointers to <i>unmutable</i> objects.
 */
//---------------------------------------------------------------------------//

class SPref
{

public:

  /// Constructor
  SPref(int r = 1) : b_refs(r) {}

  /// Return the reference count
  int refs() const
  {
    return b_refs;
  }

  /// Increment the reference count
  inline void increment()
  {
    #pragma omp critical(referencecount)
    {
      b_refs++;
    }
  }

  /// Decrement the reference count
  inline void decrement()
  {
    #pragma omp critical(referencecount)
    {
      b_refs--;
    }
  }

private:

  /// Number of references.
  int b_refs;

#ifdef DETRAN_ENABLE_BOOST

  /// Serialize
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & b_refs;
  }

#endif

};

//---------------------------------------------------------------------------//
/*!
 * \class SP
 * 
 * \brief Smart pointer implementation that does reference counting.
 *
 * The smart pointer provides a "safe" encapsulation for a standard C++
 * pointer.  Consider: A function new's an object and return the pointer to
 * that object as its return value.  Now it is the caller's responsibility to
 * free the object.  What if the caller passes the pointer to other objects
 * or functions?  What if it is not known which will be deleted first or
 * last?
 *
 * Instead the function can return a "smart pointer".  This SP class uses
 * reference counting to determine the number of current users of a pointer.
 * Each time an SP goes out of scope, the reference count is decremented.
 * When the last user of a pointer is done, the pointer is freed.
 *
 * Note: I am calling this a "smart pointer", not a "safe pointer".  There
 * are clearly ways you can hose this.  In particular, when you bind an SP<T>
 * to a T*, you yield all rights to the T*.  You'd better not squirrel the
 * bare pointer away somewhere and expect to clandestinely use it in other
 * ways or places--death will be sure to follow.  Consequently then, the
 * safest way to use this smart pointer, is to bind it to the contained
 * pointer and then always use the smart pointer.  Immediately returning the
 * smart pointer as a return value, allowing the original bare pointer to go
 * out of scope never to be seen again, is one good example of how to use
 * this.
 * 
 * One good example of bad usage is assigning the same dumb pointer to
 * multiple SPs.  Consider:
 * \code
 *     SP<Foo> f1;
 *     SP<Foo> f2;
 *     // ...
 *     Foo *f = new Foo;
 *     f1 = f;
 *     // ...
 *     f2 = f; // bad, now f1 and f2 assume they "own" f!
 * \endcode
 * Unfortunately, there is no way to check if another SP owns the dumb
 * pointer that you give to a SP.  This is simply something that needs to be
 * watched by the programmer.
 *
 * Note, SP is minimally thread-safe in that worker threads can make
 * copies of the SP (which increments the counter) and when out of scope,
 * the counter is decremented, all safely in the counter.  The objects
 * to which these SP's point are <b>not</b> thread safe, but there are
 * few, if any, cases where that behavior would be required.
 *
 */
/*!
 *  \example utilities/test/test_SP.cc
 *
 *  Test of class SP.
 */
//---------------------------------------------------------------------------//

template<class T>
class SP 
{

public:

  /// Default constructor.
  SP() : p(0), r(new SPref) { Ensure (r); Ensure (r->refs() == 1); }

  /// Explicit constructor for type T *.
  inline explicit SP(T *p_in);

  // Explicit constructor for type X *.
  template<class X>
  inline explicit SP(X *px_in);

  // Copy constructor for SP<T>.
  inline SP(const SP<T> &sp_in);

  // Copy constructor for SP<X>.
  template<class X>
  inline SP(const SP<X> &spx_in);

  /// Destructor, memory is released when count goes to zero.
  ~SP() { free(); }

  // Assignment operator for type T *.
  inline SP<T>& operator=(T *p_in);

  // Assignment operator for type X *.
  template<class X>
  inline SP<T>& operator=(X *px_in);

  // Assignment operator for type SP<T>.
  inline SP<T>& operator=(const SP<T> sp_in);

  // Assignment operator for type SP<X>.
  template<class X>
  inline SP<T>& operator=(const SP<X> spx_in);

  /// Access operator.
  T* operator->() const
  {
    Requirev(p, std::string(typeid(T).name()));
    return p;
  }

  /// Dereference operator.
  T& operator*() const
  {
    Requirev(p, std::string(typeid(T).name()));
    return *p;
  }

  /// Get the base-class pointer; better know what you are doing.
  T* bp() const { return p; }

  /// Boolean conversion operator.
  operator bool() const { return p != 0; }

  /// Operator not.
  bool operator!() const { return p == 0; }

  /// Equality operator for T*.
  bool operator==(const T *p_in) const { return p == p_in; }

  /// Inequality operator for T*.
  bool operator!=(const T *p_in) const { return p != p_in; }

  /// Equality operator for SP<T>.
  bool operator==(const SP<T> &sp_in) const { return p == sp_in.p; }

  /// Inequality operator for SP<T>.
  bool operator!=(const SP<T> &sp_in) const { return p != sp_in.p; }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Raw pointer held by smart pointer.
  T *p;

  /// Pointer to reference counter.
  SPref *r;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Free the pointer.
  inline void free();

  /// All derivatives of SP are friends.
  template<class X> friend class SP;

#ifdef DETRAN_ENABLE_BOOST

  friend class boost::serialization::access;
  template<class Archive>

  void serialize(Archive & ar, const unsigned int version)
  {
    ar & p;
    ar & r;
  }

#endif

};

} // end namespace detran_utilities

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE MEMBERS
//---------------------------------------------------------------------------//

#include "SP.i.hh"

//---------------------------------------------------------------------------//
// SP CONSTRUCTOR MACRO
//---------------------------------------------------------------------------//

#define SPCREATE(class_name, return_type, types, vals)  \
static return_type CREATE types                         \
{                                                       \
  return_type p(new class_name vals) ;                  \
  return p;                                             \
};

#endif // SP_HH_

//---------------------------------------------------------------------------//
//              end of SP.hh
//---------------------------------------------------------------------------//
