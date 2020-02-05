//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DDF.hh
 *  @brief DDF class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_DDF_HH_
#define detran_orthog_DDF_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class DDF
 *  @brief Discrete delta function basis
 *
 *  This initial implementation just piggy-backs on the
 *  default methods.  This could be improved with specialized
 *  fold and transform methods and their inverses.
 */
class DDF: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  DDF(const Parameters &p);

  /// Virtual destructor
  virtual ~DDF(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DDF_HH_

//----------------------------------------------------------------------------//
//              end of file DDF.hh
//----------------------------------------------------------------------------//
