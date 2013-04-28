//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MomentToDiscrete.i.hh
 *  @author Jeremy Roberts
 *  @date   Sep 6, 2011
 *  @brief  MomentToDiscrete inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_MOMENT_TO_DISCRETE_I_HH_
#define detran_angle_MOMENT_TO_DISCRETE_I_HH_

#include "SphericalHarmonics.hh"
#include "utilities/Constants.hh"
#include "utilities/DBC.hh"
#include "utilities/SoftEquivalence.hh"
#include <cmath>
#include <iostream>

namespace detran_angle
{

//---------------------------------------------------------------------------//
inline void MomentToDiscrete::build(SP_quadrature q)
{
  Require(q);
  Require(q->dimension() >= 1 && q->dimension() <= 3);

  // Store the angular mesh.
  d_quadrature = q;

  // Clear the existing M matrix.
  d_M.clear();
  Assert(d_M.empty());

  // Resize M.
  d_number_angles = d_quadrature->number_angles();
  Ensure(d_number_angles > 0);
  d_M.resize(d_number_angles, M_Row(d_number_moments));

  // Build the moment-to-discrete operator by looping through each octant
  // and then the angles in each octant.  The ordering is determined by the
  // order of octants in the quadrature.
  for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
  {
    for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
    {
      if (d_quadrature->dimension() == 1)
        calc_row_1d(o, a);
      else if (d_quadrature->dimension() == 2)
        calc_row_2d(o, a);
      else
        calc_row_3d(o, a);
    }
  }

}

//---------------------------------------------------------------------------//
// 3-d
inline void MomentToDiscrete::calc_row_3d(const size_t o, const size_t a)
{
  int angle = d_quadrature->index(o, a);
  // reference to current row
  M_Row &row = d_M[angle];
  // calculate the moments and add them to the row
  double mu  = d_quadrature->mu(o, a);
  double eta = d_quadrature->eta(o, a);
  double xi  = d_quadrature->xi(o, a);
  for (size_t i = 0; i < d_number_moments; ++i)
  {
    size_t l = d_indexer->l(i);
    int    m = d_indexer->m(i);
    double norm = (2.0 * l + 1.0) * detran_utilities::inv_four_pi;
    row[i] = norm * SphericalHarmonics::Y_lm(l, m, mu, eta, xi);
  }
}

//---------------------------------------------------------------------------//
// 2-d
inline void MomentToDiscrete::calc_row_2d(const size_t o, const size_t a)
{
  int angle = d_quadrature->index(o, a);
  // reference to current row
  M_Row &row = d_M[angle];
  // compute the direction cosine w/r to polar axis
  double mu  = d_quadrature->mu(o, a);
  double eta = d_quadrature->eta(o, a);
  double xi  = std::sqrt(1.0 - mu * mu - eta * eta);
  for (size_t i = 0; i < d_number_moments; ++i)
  {
    size_t l = d_indexer->l(i);
    int    m = d_indexer->m(i);
    double norm = (2.0 * l + 1.0) * detran_utilities::inv_four_pi;
    row[i] = norm * SphericalHarmonics::Y_lm(l, m, mu, eta, xi);
  }
}

//---------------------------------------------------------------------------//
// 1-d
inline void MomentToDiscrete::calc_row_1d(const size_t o, const size_t a)
{
  int angle = d_quadrature->index(o, a);
  // reference to current row
  M_Row &row = d_M[angle];
  // calculate the moments and add them to the row
  for (size_t l = 0; l < d_number_moments; ++l)
  {
    double norm = (2.0 * l + 1.0) * 0.5;
    row[l] = norm * SphericalHarmonics::Y_lm(l, d_quadrature->mu(o, a));
  }
}

//-----------------------------------------------------------------------------
// Operator(cardinal_angle_index, cardinal_moment_index)
inline const double& MomentToDiscrete::
operator()(const size_t angle, const size_t moment) const
{
  Require(angle < d_number_angles);
  Require(moment < d_number_moments);
  return d_M[angle][moment];
}

//-----------------------------------------------------------------------------
// Operator(cardinal_angle_index, legendre_degree, legendre_order)
inline const double& MomentToDiscrete::
operator()(const size_t angle, const size_t l, const size_t m) const
{
  // Moment cardinal index
  size_t moment = d_indexer->index(l, m);
  return (*this)(angle, moment);
}

//-----------------------------------------------------------------------------
// Operator(octant_index, angle_in_octant, legendre_degree, legendre_order)
inline const double& MomentToDiscrete::
operator()(const size_t o, const size_t a,
           const size_t l, const size_t m) const
{
  // Angle cardinal index
  size_t angle = d_quadrature->index(o, a);
  // Moment cardinal index
  size_t moment = d_indexer->index(l, m);
  return (*this)(angle, moment);
}

} // end namespace detran_angle

#endif /* detran_angle_MOMENT_TO_DISCRETE_I_HH_ */
