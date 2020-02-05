//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.t.hh
 *  @brief ProductQuadrature member definitions
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_PRODUCTQUADRATURE_T_HH_
#define detran_angle_PRODUCTQUADRATURE_T_HH_

#include "ProductQuadratureAP.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
template <class A, class P>
ProductQuadratureAP<A, P>::ProductQuadratureAP(const size_t dim,
                                               const size_t na,
                                               const size_t np,
                                               const bool   normalize)
  : ProductQuadrature(dim, na, np,
                      "product-" + A::name() + "-" + P::name(), normalize)
{
  // build azimuth
  A Q_A(d_number_azimuth_octant, normalize);
  d_phi = Q_A.phi();
  d_cos_phi = Q_A.cos_phi();
  d_sin_phi = Q_A.sin_phi();
  d_azimuth_weight = Q_A.weights();

  // build polar
  P Q_P(d_number_polar_octant, normalize);
  d_cos_theta = Q_P.cos_theta();
  d_sin_theta = Q_P.sin_theta();
  d_polar_weight = Q_P.weights();

  build();
}

} // end namespace detran_angle

#endif /* detran_angle_PRODUCTQUADRATURE_T_HH_ */

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.t.hh
//----------------------------------------------------------------------------//
