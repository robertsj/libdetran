//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CLP.hh
 *  @brief CLP class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_CLP_HH_
#define detran_orthog_CLP_HH_

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class CLP
 *  @brief Continuous Legendre polynomial basis
 */
class ORTHOG_EXPORT CLP: public ContinuousOrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  CLP(const Parameters &p);

  /// Virtual destructor
  virtual ~CLP(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_CLP_HH_

//----------------------------------------------------------------------------//
//              end of file CLP.hh
//----------------------------------------------------------------------------//
