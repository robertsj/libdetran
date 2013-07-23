//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DCT.hh
 *  @brief DCT class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_DCT_HH_
#define detran_orthog_DCT_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class DCT
 *  @brief Discrete cosine transform
 */
class ORTHOG_EXPORT DCT: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  DCT(const Parameters &p);

  /// Virtual destructor
  virtual ~DCT(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DCT_HH_

//----------------------------------------------------------------------------//
//              end of file DCT.hh
//----------------------------------------------------------------------------//
