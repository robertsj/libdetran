//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PCShell.hh
 *  @author robertsj
 *  @date   Sep 20, 2012
 *  @brief  PCShell class definition.
 */
//---------------------------------------------------------------------------//

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

class PCShell: public Preconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Preconditioner                   Base;
  typedef typename Base::SP_preconditioner SP_preconditioner;
  typedef typename Vector::SP_vector       SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Construct a shell preconditioner
  PCShell(std::string name = "PCShell");

  /// Virtual destructor
  virtual ~PCShell(){};

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /* ... */

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = b
  virtual void apply(Vector &b, Vector &x) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  using Base::d_name;

};

} // end namespace callow

#endif /* callow_PCSHELL_HH_ */
