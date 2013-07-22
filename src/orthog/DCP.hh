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
 */
class ORTHOG_EXPORT DCP: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *   @brief Constructor.
   *   @param   order   Order of the basis
   *   @param   size    Size of the basis vectors
   */
  DCP(const size_t order, const size_t size);

  /// Virtual destructor
  virtual ~DCP(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DCP_HH_

//----------------------------------------------------------------------------//
//              end of file DCP.hh
//----------------------------------------------------------------------------//
