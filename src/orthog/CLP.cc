//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CLP.cc
 *  @brief  CLP
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#include "CLP.hh"
#ifdef DETRAN_ENABLE_BOOST
#include <boost/math/special_functions/legendre.hpp>
#endif

namespace detran_orthog
{

//---------------------------------------------------------------------------//
CLP::CLP(const size_t order,
         const vec_dbl &x,
         const vec_dbl &dx)
  : ContinuousOrthogonalBasis(order, x, dx)
{
#ifndef DETRAN_ENABLE_BOOST
  THROW("CLP needs boost to be enabled.");
#else

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_order + 1, d_size, 0.0);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // The weights are just the dx's.
  for (size_t i = 0; i < d_w->size(); ++i)
    (*d_w)[i] = d_dx[i];

  // Compute left and right edges.  The domain must be translated
  // and scaled to be in [-1, 1]
  double a = d_x[0] - 0.5*d_dx[0];
  double b = d_x[d_size - 1] + 0.5*d_dx[d_size - 1];
  double L = b - a;

  // Build the basis
  for (size_t l = 0; l <= d_order; ++l)
  {
    for (size_t i = 0; i < d_size; ++i)
    {
      // Scale the abscissa
      double z =  2.0*(d_x[i] - a)/L - 1;
      (*d_basis)(l, i) = boost::math::legendre_p(l, z);
    }
    // Inverse of normalization coefficient.
    (*d_a)[l] = (2.0 * l + 1.0) / 2.0;
  }

#endif
}

} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file CLP.cc
//---------------------------------------------------------------------------//
