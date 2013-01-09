//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DCT.cc
 *  @brief  DCT
 *  @author Jeremy Roberts
 *  @date   Jan 9, 2013
 */
//---------------------------------------------------------------------------//

#include "DCT.hh"
#include "utilities/Constants.hh"
#include <cmath>

namespace detran_orthog
{

//---------------------------------------------------------------------------//
DCT::DCT(const size_t order, const size_t size)
  : OrthogonalBasis(order, size)
{
  using detran_utilities::pi;

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_order + 1, d_size, 1.0);

  for (int i = 0; i <= d_order; ++i)
    for (int j = 0; j < d_size; ++j)
      (*d_basis)(i, j) = std::cos((pi/d_size) * (j + 0.5) * i);
  compute_a();

}


} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file DCT.cc
//---------------------------------------------------------------------------//
