//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DLP.hh
 *  @brief  DLP
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_DLP_HH_
#define detran_orthog_DLP_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class DLP
 *  @brief Discrete Legendre polynomial basis
 */
class DLP: public OrthogonalBasis
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *   @brief Constructor.
   *   @param   order       Order of the basis
   *   @param   size        Size of the basis vectors
   *   @param   orthonormal Flag to indicate the basis should be orthonormal
   */
  DLP(const size_t order,
      const size_t size,
      const bool orthonormal = false);

  /// Virtual destructor
  virtual ~DLP(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DLP_HH_

//---------------------------------------------------------------------------//
//              end of file DLP.hh
//---------------------------------------------------------------------------//
