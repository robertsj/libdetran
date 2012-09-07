//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MomentToDiscrete.i.hh
 * \author Jeremy Roberts
 * \date   Sep 6, 2011
 * \brief  MomentToDiscrete inline member definitions.
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef MOMENT_TO_DISCRETE_I_HH_
#define MOMENT_TO_DISCRETE_I_HH_

// Detran
#include "SphericalHarmonics.hh"

// Utilities
#include "Constants.hh"
#include "DBC.hh"
#include "SoftEquivalence.hh"

// System
#include <cmath>
#include <iostream>

namespace detran_angle
{

// Constructor
MomentToDiscrete::MomentToDiscrete(const size_t legendre_order)
  : d_legendre_order(legendre_order)
  , d_number_moments(1) // \todo save for later use
{
  Require(d_legendre_order >= 0);
  Require(d_number_moments > 0);
}

// Build the operator.
void MomentToDiscrete::build(SP_quadrature q)
{
  Require(q);
  Require(q->dimension() >= 1 and q->dimension() <= 3);

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
  for (int o = 0; o < d_quadrature->number_octants(); o++)
  {
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
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

// 3-d
inline void MomentToDiscrete::calc_row_3d(const size_t o, const size_t a)
{
  Assert(d_quadrature);
  int angle = d_quadrature->index(o, a);

  // reference to current row
  M_Row &row = d_M[angle];

  // calculate the moments and add them to the row
  for (int l = 0; l <= d_legendre_order; l++)
  {
    double norm = (2.0 * l + 1.0) * detran_utilities::inv_four_pi;
    ;
    for (int m = -l; m <= l; m++)
    {
      // get cardinal moment index
      int i = 0; //Moments<D>::index(l,m);
      Assert(i < d_number_moments);
      // add harmonic
      row[i] = norm
          * SphericalHarmonics::Y_lm(l, m,
                                     d_quadrature->mu(o, a),
                                     d_quadrature->eta(o, a),
                                     d_quadrature->xi(o, a));
    }
  }
}

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

  // loop through l>0 moments and add
  for (int l = 0; l <= d_legendre_order; l++)
  {
    double norm = (2.0 * l + 1.0) * detran_utilities::inv_four_pi;
    for (int m = -l; m <= l; m += 2)
    {
      int i = 0; // Moments<_2D>::index(l, m);
      row[i] = norm * SphericalHarmonics::Y_lm(l, m, mu, eta, xi);
    }
  }
}

// 1-d
inline void MomentToDiscrete::calc_row_1d(const size_t o, const size_t a)
{
  int angle = d_quadrature->index(o, a);

  // reference to current row
  M_Row &row = d_M[angle];

  // calculate the moments and add them to the row
  for (int l = 0; l <= d_legendre_order; l++)
  {
    double norm = (2.0 * l + 1.0) * 0.5;
    int i = 0; //Moments<_1D>::index(l, 0);
    row[i] = norm * SphericalHarmonics::Y_lm(l, d_quadrature->mu(o, a));
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
  int index = 0; //Moments<D>::index(l, m);
  return (*this)(angle, index);
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
  size_t moment = 0; //Moments<D>::index(l, m);
  return (*this)(angle, moment);
}

} // end namespace detran_angle

#endif /* MOMENT_TO_DISCRETE_I_HH_ */
