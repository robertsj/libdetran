//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DCT.cc
 *  @brief DCT member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "DCT.hh"
#include "utilities/Constants.hh"
#include <cmath>

namespace detran_orthog
{

//----------------------------------------------------------------------------//
DCT::DCT(const Parameters &p)
  : OrthogonalBasis(p)
{
  using detran_utilities::pi;

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_size, d_order + 1, 1.0);

  for (int i = 0; i <= d_order; ++i)
    for (int j = 0; j < d_size; ++j)
      (*d_basis)(j, i) = std::cos((pi/d_size) * (j + 0.5) * i);

  compute_a();
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file DCT.cc
//----------------------------------------------------------------------------//
