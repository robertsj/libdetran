//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobi01.hh
 *  @brief Jacobi01 class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_Jacobi01_HH_
#define detran_orthog_Jacobi01_HH_

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class Jacobi01
 *  @brief Continuous Jacobi_01 polynomial basis
 *
 *  Implements the Jacobi polynomial with alpha = 0, beta = 1.
 */
class ORTHOG_EXPORT Jacobi01: public ContinuousOrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  Jacobi01(const Parameters &p);

  /// Virtual destructor
  virtual ~Jacobi01(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_Jacobi01_HH_

//----------------------------------------------------------------------------//
//              end of file Jacobi01.hh
//----------------------------------------------------------------------------//
