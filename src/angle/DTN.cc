//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DTN.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  DTN class definition.
 */
//---------------------------------------------------------------------------//

#include "DTN.hh"
#include "GenerateGaussChebyshev.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
DTN::DTN(size_t order)
  : Quadrature(1, order, "DTN")
{

  // generate the parameters on-the-fly.
  vec_dbl tmp_mu(d_number_angles_octant, 0.0);
  vec_dbl tmp_wt(d_number_angles_octant, 0.0);
  generate_gc_parameters(d_number_angles_octant, tmp_mu, tmp_wt);

  for (int i = 0; i < d_number_angles_octant; i++)
  {
    // Shift the mu to [0, 1] from [-1, 1]
    d_mu[i]     = 0.5*tmp_mu[i] + 0.5;
    d_weight[i] = tmp_wt[i];
  }

}

} // end namespace detran_angle


