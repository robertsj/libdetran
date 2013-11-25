//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DCP.cc
 *  @brief DCP member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "DCP.hh"
#include <cmath>

namespace detran_orthog
{

//----------------------------------------------------------------------------//
DCP::DCP(const Parameters &p)
  : OrthogonalBasis(p)
{
  double N = (double) d_size;

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_size, d_order + 1, 1.0);

  // Zeroth order is just a constant
  for (int j = 0; j < d_size; ++j)
    (*d_basis)(j, 0) = 1.0;

  // Do others, using symmetry
  size_t mid = d_size % 2 ? mid = (d_size + 1)/2 : mid = d_size / 2;
  for (int i = 1; i <= d_order; ++i)
  {
    double g = -std::sqrt((N-i)/(N+i)) *
                std::sqrt((2.0*i + 1.0) / (2.0*i - 1.0));
    (*d_basis)(0, i) = g * (*d_basis)(0, i - 1);
    g = (1.0 + (i * (i + 1.0) ) / (1.0 - N) );
    (*d_basis)(1, i)  = g *(*d_basis)(0, i);
    if (d_size > 2)
    {
      for (int j = 2; j < mid; ++j)
      {
        double g1 = (-(i*(i+1.0))-(2.0*j - 1.0)*(j-N-1.0)-j) /
                    (j*(N-j));
        double g2 = (j-1.0)*(j-N-1.0) / (j*(N-j));
        (*d_basis)(j, i) = g1*(*d_basis)(j-1, i) + g2*(*d_basis)(j-2, i);
      }
    }
    // Apply symmetry for other elements
    for (int j = d_size-1; j >= mid; --j)
      (*d_basis)(j, i) = (1 - 2 * (i % 2) ) * (*d_basis)(d_size-j-1, i);
  }
  compute_a();

}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file DCP.cc
//----------------------------------------------------------------------------//
