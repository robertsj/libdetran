//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ChebyshevU.cc
 *  @brief ChebyshevU
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ChebyshevU.hh"

namespace detran_orthog
{

//----------------------------------------------------------------------------//
ChebyshevU::ChebyshevU(const Parameters &p)
  : ContinuousOrthogonalBasis(p)
{
  size_t even = 1;
  if (d_even_only) even = 2;

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_size, d_order + 1, 0.0);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // The weights are the quadrature weights times the actual weight
  double L = d_upper_bound - d_lower_bound;
  for (size_t i = 0; i < d_w->size(); ++i)
  {
    d_x[i] = 2.0 * (d_x[i] - d_lower_bound) / L - 1.0;
    (*d_w)[i] = d_qw[i] * std::sqrt(1.0 - d_x[i] * d_x[i]);
  }

  // Build the basis
  for (size_t l = 0; l <= d_order; ++l)
    for (size_t i = 0; i < d_size; ++i)
      (*d_basis)(i, l) = cheby_u(even * l, d_x[i]);


  //d_orthonormal = true;
  compute_a();
}

//----------------------------------------------------------------------------//
double ChebyshevU::cheby_u(const size_t l, const double x)
{
  using std::pow;
  double v = 0;
  if (l == 0)
    v = 1.0;
  else if (l == 1)
    v = 2.0 * x;
  else if (l == 2)
    v = -1.0 + 4.0 * x * x;
  else if (l == 3)
    v = 0.8e1 * pow(x, 0.3e1) - 0.4e1 * x;
  else if (l == 4)
    v = 0.1e1 + 0.16e2 * pow(x, 0.4e1) - 0.12e2 * x * x;
  else if (l == 5)
    v = 0.32e2 * pow(x, 0.5e1) - 0.32e2 * pow(x, 0.3e1) + 0.6e1 * x;
  else if (l == 6)
    v = -0.1e1 + 0.64e2 * pow(x, 0.6e1) - 0.80e2 * pow(x, 0.4e1)
        + 0.24e2 * x * x;
  else if (l == 7)
    v = 0.128e3 * pow(x, 0.7e1) - 0.192e3 * pow(x, 0.5e1)
        + 0.80e2 * pow(x, 0.3e1) - 0.8e1 * x;
  else if (l == 8)
    v = 0.1e1 + 0.256e3 * pow(x, 0.8e1) - 0.448e3 * pow(x, 0.6e1)
        + 0.240e3 * pow(x, 0.4e1) - 0.40e2 * x * x;
  else
    v = 2.0 * cheby_u(l - 1, x) - cheby_u(l - 2, x);
  return v;
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file ChebyshevU.cc
//----------------------------------------------------------------------------//
