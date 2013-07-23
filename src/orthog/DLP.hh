//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DLP.hh
 *  @brief DLP class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_DLP_HH_
#define detran_orthog_DLP_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class DLP
 *  @brief Discrete Legendre polynomial basis
 */
class ORTHOG_EXPORT DLP: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  DLP(const Parameters &p);

  /// Virtual destructor
  virtual ~DLP(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DLP_HH_

//----------------------------------------------------------------------------//
//              end of file DLP.hh
//----------------------------------------------------------------------------//
