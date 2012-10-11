//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussChebyshev.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  GaussChebyshev class definition.
 */
//---------------------------------------------------------------------------//

#include "GaussChebyshev.hh"
#include "GenerateGaussChebyshev.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
GaussChebyshev::GaussChebyshev(size_t order)
  : Quadrature(1, order, "GaussChebyshev")
{
  // Preconditions
  Require(order % 2 == 0);

  // generate the parameters on-the-fly.
  vec_dbl tmp_mu(order, 0.0);
  vec_dbl tmp_wt(order, 0.0);
  generate_gc_parameters(order, tmp_mu, tmp_wt, false);

  for (int i = 0; i < d_number_angles_octant; i++)
  {
    // Shift the mu to [0, 1] from [-1, 1]
    d_mu[i]     = 0.5*tmp_mu[i] + 0.5;
    d_weight[i] = tmp_wt[i];
  }

}

} // end namespace detran_angle



