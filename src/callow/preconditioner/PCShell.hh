//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCShell.hh
 *  @brief PCShell class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_PCSHELL_HH_
#define callow_PCSHELL_HH_

#include "Preconditioner.hh"

namespace callow
{

/**
 *  @class PCShell
 *  @brief Applies a shell preconditioner
 *
 *  A shell preconditioner allows the user to define their own
 *  preconditioning processes that are potentially matrix free.
 */

class CALLOW_EXPORT PCShell: public Preconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Preconditioner                Base;
  typedef Base::SP_preconditioner       SP_preconditioner;
  typedef Vector::SP_vector             SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Construct a shell preconditioner
  PCShell(std::string name = "PCShell", void* context = NULL);

  /// Virtual destructor
  virtual ~PCShell(){};

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /* ... */

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// Solve Px = b
  virtual void apply(Vector &b, Vector &x) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// User context
  void* d_context;

};

} // end namespace callow

#endif /* callow_PCSHELL_HH_ */

//----------------------------------------------------------------------------//
//              end of PCShell.hh
//----------------------------------------------------------------------------//
