//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  UserBasis.hh
 *  @brief UserBasis class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_UserBasis_HH_
#define detran_orthog_UserBasis_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class UserBasis
 *  @brief User defined discrete functional basis
 */
class ORTHOG_EXPORT UserBasis: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  UserBasis(const Parameters &p);

  /// Virtual destructor
  virtual ~UserBasis(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_UserBasis_HH_

//----------------------------------------------------------------------------//
//              end of file UserBasis.hh
//----------------------------------------------------------------------------//
