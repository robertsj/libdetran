//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DCP.hh
 *  @brief DCP class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_DCP_HH_
#define detran_orthog_DCP_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class DCP
 *  @brief Discrete Chebyshev polynomial basis
 *
 *  Note, it turns out these are equivalent to @ref DLP.  Basically, given
 *  a flat and linear starting basis, there can only be one next step
 *  in the recursion (outside of normalization).  However, this
 *  implementation is slightly more stable numerically.
 */
class ORTHOG_EXPORT DCP: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  DCP(const Parameters &p);

  /// Virtual destructor
  virtual ~DCP(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DCP_HH_

//----------------------------------------------------------------------------//
//              end of file DCP.hh
//----------------------------------------------------------------------------//
