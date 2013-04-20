//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DLP.cc
 *  @brief  DLP
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#include "DLP.hh"

namespace detran_orthog
{

//---------------------------------------------------------------------------//
DLP::DLP(const size_t   order,
         const size_t   size,
         const bool     orthonormal)
  : OrthogonalBasis(order, size, orthonormal)
{

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_order + 1, d_size, 1.0);

  // Zeroth built by initialization. Build first order if needed.
  if (d_order > 0)
  {
    for (int j = 0; j < d_size; ++j)
      (*d_basis)(1, j) = (d_size - 1 - (2.0 * j))/(d_size - 1);
    Vector row(d_size, &(*d_basis)(1, 0));

    // Higher orders
    for (int i = 2; i <= d_order; ++i)
    {
      for (int j = 0; j < d_size; ++j)
      {
        int C0 = (i - 1) * (d_size - 1 + i);
        int C1 = (2 * i - 1) * (d_size - 1 - 2 * j);
        int C2 = i * (d_size - i);
        (*d_basis)(i, j) = (C1 * (*d_basis)(i-1, j) -
                            C0 * (*d_basis)(i-2, j) ) / C2;
      }
    }
  }
  compute_a();
}


} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file DLP.cc
//---------------------------------------------------------------------------//
