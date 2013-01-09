//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DCT.hh
 *  @brief  DCT
 *  @author Jeremy Roberts
 *  @date   Jan 9, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_DCT_HH_
#define detran_orthog_DCT_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{


/**
 *  @class DCT
 *  @brief Discrete cosine transform
 */
class DCT: public OrthogonalBasis
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
  DCT(const size_t order, const size_t size);

  /// Virtual destructor
  virtual ~DCT(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DCT_HH_

//---------------------------------------------------------------------------//
//              end of file DCT.hh
//---------------------------------------------------------------------------//
