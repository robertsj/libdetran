//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Jacobi01.cc
 *  @brief  Jacobi01
 *  @author Jeremy Roberts
 *  @date   Jan 21, 2013
 */
//---------------------------------------------------------------------------//

#include "Jacobi01.hh"

namespace detran_orthog
{

double jacobi_coefs [8][8] =
 {
  {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {-0.5, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  {-0.5, -1.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0},
  {0.375, -1.875, -1.875, 4.375, 0.0, 0.0, 0.0, 0.0},
  {0.375, 1.5, -5.25, -3.5, 7.875, 0.0, 0.0, 0.0},
  {-0.3125, 2.1875, 4.375, -13.1250, -6.5625, 14.4375, 0.0, 0.0},
  {-0.3125, -1.875, 8.4375, 11.25, -30.9375, -12.3750, 26.8125, 0.0},
  {0.2734375, -2.4609375, -7.3828125, 27.07031250, 27.0703125, -70.3828125, -23.4609375, 50.27343750}
 };



//---------------------------------------------------------------------------//
Jacobi01::Jacobi01(const size_t order,
                   const vec_dbl &x,
                   const vec_dbl &qw,
                   const double x_0,
                   const double x_1)
  : ContinuousOrthogonalBasis(order, x, qw)
{
  Insist(order < 8, "Maximum of 8th order for Jacobi01.");

  d_orthonormal = true;

  // Allocate the basis matrix
  d_basis = new callow::MatrixDense(d_order + 1, d_size, 0.0);

  // Allocate the normalization array
  d_a = Vector::Create(d_order + 1, 0.0);

  // The weights are the quadrature weights times the actual weight
  double L = x_1 - x_0;
  for (size_t i = 0; i < d_w->size(); ++i)
  {
    // Scale the abscissa 2*[0 1] - 1    2 * (x(:)+x_0)/L
    d_x[i] =  2.0*(d_x[i] - x_0)/L - 1.0;
    (*d_w)[i] = d_qw[i] * 0.5 * (d_x[i] + 1.0);
  }

  // Build the basis
  for (size_t l = 0; l <= d_order; ++l)
  {
    for (size_t i = 0; i < d_size; ++i)
    {
      for (size_t j = 0; j <= l; ++j)
      {
        (*d_basis)(l, i) += jacobi_coefs[l][j] * std::pow(d_x[i], j) *
                            std::sqrt(l+1);
      }
    }
  }
  compute_a();
}

} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file Jacobi01.cc
//---------------------------------------------------------------------------//
