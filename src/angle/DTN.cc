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
DTN::DTN(const size_t number_polar,
         const bool   normalize)
  : Quadrature(1, 2*number_polar, "DTN")
{

  // generate the parameters on-the-fly.
  vec_dbl tmp_mu(d_number_angles_octant, 0.0);
  vec_dbl tmp_wt(d_number_angles_octant, 0.0);
  generate_gc_parameters(d_number_angles_octant, tmp_mu, tmp_wt, normalize);

  for (size_t i = 0; i < d_number_angles_octant; ++i)
  {
    // Shift the mu to [0, 1] from [-1, 1] and halve the weight
    d_mu[i]     = 0.5*tmp_mu[i] + 0.5;
    d_weight[i] = 0.5*tmp_wt[i];
  }

}

//---------------------------------------------------------------------------//
DTN::SP_quadrature
DTN::Create(const size_t number_polar_octant,
            const bool   normalize)
{
  SP_quadrature p(new DTN(number_polar_octant, normalize));
  return p;
}

} // end namespace detran_angle


