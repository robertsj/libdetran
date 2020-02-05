//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DLP.cc
 *  @brief DLP member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "DLP.hh"

namespace detran_orthog
{

//----------------------------------------------------------------------------//
DLP::DLP(const Parameters &p)
  : OrthogonalBasis(p)
{
  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_size, d_order + 1, 1.0);

  // Zeroth built by initialization. Build first order if needed.
  if (d_order > 0)
  {
    for (int j = 0; j < d_size; ++j)
      (*d_basis)(j, 1) = (d_size - 1 - (2.0 * j))/(d_size - 1);

    // Higher orders
    for (int i = 2; i <= d_order; ++i)
    {
      for (int j = 0; j < d_size; ++j)
      {
        int C0 = (i - 1) * (d_size - 1 + i);
        int C1 = (2 * i - 1) * (d_size - 1 - 2 * j);
        int C2 = i * (d_size - i);
        (*d_basis)(j, i) = (C1 * (*d_basis)(j, i-1) -
                            C0 * (*d_basis)(j, i-2) ) / C2;
      }
    }
  }

  compute_a();
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file DLP.cc
//----------------------------------------------------------------------------//
