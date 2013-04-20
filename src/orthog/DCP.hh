//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DCP.hh
 *  @brief  DCP
 *  @author Jeremy Roberts
 *  @date   Jan 9, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_DCP_HH_
#define detran_orthog_DCP_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{


/**
 *  @class DCP
 *  @brief Discrete Chebyshev polynomial basis
 */
class DCP: public OrthogonalBasis
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

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

//---------------------------------------------------------------------------//
//              end of file DCP.hh
//---------------------------------------------------------------------------//
