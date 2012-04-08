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

namespace detran
{

  // Constructor.
  template <class D>
  MomentToDiscrete<D>::MomentToDiscrete(const size_type legendre_order)
    :  d_legendre_order(legendre_order)
    ,  d_number_moments(1)             // \todo save for later use
  {
    Require(d_legendre_order >= 0);
    Require(d_number_moments > 0);
  }

  // Build the operator.
  template <class D>
  void MomentToDiscrete<D>::build(SP_quadrature quadrature)
  {
    Require(quadrature);

    // Store the angular mesh.
    d_quadrature = quadrature;

    // Clear the existing M matrix.
    d_M.clear();
    Assert(d_M.empty());

    // Resize M.
    d_number_angles = quadrature->number_angles();
    Ensure(d_number_angles > 0);
    d_M.resize(d_number_angles, M_Row(d_number_moments));
    Assert(d_M.size() == quadrature->number_angles());

    // Build the moment-to-discrete operator by looping through each octant
    // and then the angles in each octant.  The ordering is determined by the
    // order of octants in the quadrature.
    for (int o = 0; o < d_quadrature->number_octants(); o++)
    {
      for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
      {
        calc_row(o, a);
      }
    }

  }

  // 3-d
  template <class D>
  void MomentToDiscrete<D>::calc_row(int o, int a)
  {
    Assert(d_quadrature);
    int angle = d_quadrature->index(o, a);

    // reference to current row
    M_Row &row = d_M[angle];

    // calculate the moments and add them to the row
    for (int l = 0; l <= d_legendre_order; l++)
    {
        double norm = (2.0 * l + 1.0) * detran_utils::inv_four_pi;;
        for (int m = -l; m <= l; m++)
        {
            // get cardinal moment index
            int i = 0; //Moments<D>::index(l,m);
            Assert(i < d_number_moments);
            // add harmonic
            row[i] = norm * SphericalHarmonics::
                     Y_lm(l, m,
                          d_quadrature->mu(o, a),
                          d_quadrature->eta(o, a),
                          d_quadrature->xi(o, a));
        }
    }
  }

  // 2-d
  template <>
  void MomentToDiscrete<_2D>::calc_row(const int o, const int a)
  {
    int angle = d_quadrature->index(o, a);

    // reference to current row
    M_Row &row = d_M[angle];

    // compute the direction cosine w/r to polar axis
    double mu  = d_quadrature->mu(o, a);
    double eta = d_quadrature->eta(o, a);
    double xi  = sqrt(1.0 - mu*mu - eta*eta);

    // loop through l>0 moments and add
    for (int l = 0; l <= d_legendre_order; l++)
    {
      double norm = (2.0 * l + 1.0) * detran_utils::inv_four_pi;
      for (int m = -l; m <= l; m += 2)
      {
        int i = 0; // Moments<_2D>::index(l, m);
        row[i] = norm * SphericalHarmonics::
            Y_lm(l, m, mu, eta, xi);
      }
    }
  }

  // 1-d
  template <>
  void MomentToDiscrete<_1D>::calc_row(const int a, const int o)
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
template <class D>
inline const double& MomentToDiscrete<D>::
operator()(const int angle, const int moment) const
{
  Require(angle >= 0);
  Require(angle < d_number_angles);
  Require (moment >= 0);
  Require(moment < d_number_moments);
  return d_M[angle][moment];
}

//-----------------------------------------------------------------------------
// Operator(cardinal_angle_index, legendre_degree, legendre_order)
template <class D>
inline const double& MomentToDiscrete<D>::
operator()(const int angle, const int l, const int m) const
{
  // Moment cardinal index
  int index = 0; //Moments<D>::index(l, m);
  return (*this)(angle, index);
}

//-----------------------------------------------------------------------------
// Operator(octant_index, angle_in_octant, legendre_degree, legendre_order)
template <class D>
inline const double& MomentToDiscrete<D>::
operator()(const int o, const int n,
           const int l, const int m) const
{
  // Angle cardinal index
  int angle = d_quadrature->index(o, n);
  // Moment cardinal index
  int moment = 0; //Moments<D>::index(l, m);
  return (*this)(angle, moment);
}

// explicit instantiations
template class MomentToDiscrete<_1D>;
template class MomentToDiscrete<_2D>;
template class MomentToDiscrete<_3D>;

} // end namespace detran

#endif /* MOMENT_TO_DISCRETE_I_HH_ */
