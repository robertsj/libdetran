//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CLP.cc
 *  @brief CLP member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "CLP.hh"
#include "utilities/DBC.hh"

#ifdef DETRAN_ENABLE_BOOST
#include <boost/math/special_functions/legendre.hpp>
#define LP(l, x) boost::math::legendre_p(l, d_x[i]);
#else
double LP(int l, double x)
{
  if (l == 0)
    return 1.0;
  if (l == 1)
    return x;
  // (n+1)P_n+1(x) = (2n+1)xP_n(x) - nP_n-1(x)
  return ((2*(l-1)+1)*x*LP(l-1, x) - (l-1)*LP(l-2, x))/(l);
}
#endif

namespace detran_orthog
{

//----------------------------------------------------------------------------//
CLP::CLP(const Parameters &p)
  : ContinuousOrthogonalBasis(p)
{
  // Allocate the basis matrix
  d_basis  = std::make_shared<callow::MatrixDense>(d_size, d_order + 1, 0.0);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // The weights are just the qw's.
  double L = d_upper_bound - d_lower_bound;
  for (size_t i = 0; i < d_w->size(); ++i)
  {
    d_x[i] =  2.0*(d_x[i] - d_lower_bound)/L - 1.0;
    (*d_w)[i] = d_qw[i];
  }

  // Build the basis
  for (size_t l = 0; l <= d_order; ++l)
  {
    for (size_t i = 0; i < d_size; ++i)
    {
      (*d_basis)(i, l) = LP(l, d_x[i]);
    }
    // Inverse of normalization coefficient.
    (*d_a)[l] = (2.0 * l + 1.0) / 2.0;
  }

}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file CLP.cc
//----------------------------------------------------------------------------//
