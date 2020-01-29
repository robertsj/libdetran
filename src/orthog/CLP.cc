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
#endif

namespace detran_orthog
{

//----------------------------------------------------------------------------//
CLP::CLP(const Parameters &p)
  : ContinuousOrthogonalBasis(p)
{
#ifndef DETRAN_ENABLE_BOOST
  THROW("CLP needs boost to be enabled.");
#else

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_size, d_order + 1, 0.0);

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
      (*d_basis)(i, l) = boost::math::legendre_p(l, d_x[i]);
    }
    // Inverse of normalization coefficient.
    (*d_a)[l] = (2.0 * l + 1.0) / 2.0;
  }
  //d_orthonormal = true;
  //compute_a();
#endif

}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file CLP.cc
//----------------------------------------------------------------------------//
