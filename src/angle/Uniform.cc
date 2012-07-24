//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Uniform.cc
 * \brief  Uniform 
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "Uniform.hh"

// Utilities
#include "SoftEquivalence.hh"

// System
#include <cmath>
#include <iostream>


namespace detran
{

Uniform::Uniform(int dim,
                 int num_azimuths_octant,
                 int num_space,
                 int num_polar,
                 std::string polar)
  : QuadratureMOC(dim,
                  num_azimuths_octant,
                  num_polar,
                  "UNIFORM",
                  polar)
{

  using std::cout;
  using std::endl;

  bool db = false;

  // Define evenly-spaced azimuthal angles over [0, pi/2].  These
  // are subject change.
  vec_dbl phi_quadrant(d_number_azimuths_octant, 0.0);
  double delta = 0.5 * pi / d_number_azimuths_octant;
  phi_quadrant[0] = delta/2;
  for (int a = 1; a < d_number_azimuths_octant; a++)
  {
    phi_quadrant[a] = phi_quadrant[a-1] + delta;
  }
  Ensure(soft_equiv(phi_quadrant[d_number_azimuths_octant-1], pi/2-delta/2));


  vec_int num_x(d_number_azimuths_octant, 0);
  vec_int num_y(d_number_azimuths_octant, 0);
  vec_dbl phi(d_number_azimuths_octant, 0.0);

  for (int a = 0; a < d_number_azimuths_octant; a++)
  {

    if (db) cout << " a = " << a << endl;

    double tan_phi = tan(phi_quadrant[a]);

    if (db) cout << " tan_phi = " << tan_phi << " phi = " << phi_quadrant[a] << endl;

    // This keeps the adjusted angles symmetric about pi/4.
    if (a < d_number_azimuths_octant / 2)
    {
        num_x[a] = std::abs(std::ceil(
          double(num_space) * tan_phi / (tan_phi + 1.0)));
    }
    else
    {
        num_x[a] = std::abs(std::floor(
          double(num_space) * tan_phi / (tan_phi + 1.0)));
    }
    num_y[a] = num_space - num_x[a];

    Assert(num_x[a] > 0);
    Assert(num_y[a] > 0);

    // cout << " num_xy[a]  = " << num_x[a] << " " << num_y[a]  << endl;

    // Actual azimuth.
    phi[a] = atan(double(num_x[a]) / double(num_y[a]));

    if (db) cout << " phi[a]  = " << phi[a]  << endl;

    // \todo Weight adjustment
    d_azimuth_weight[a] = delta;

    if (db) cout << " azimuthal_weight[a]  = " << d_azimuth_weight[a]  << endl;

    // Loop over polar angles.  Polar is on the inner, since that's
    // where it is in the sweeps.
    for (int p = 0; p < d_number_polar; p++)
    {

      // Cardinal Index
      int angle = a * d_number_polar + p;

      // Fill the angles.  Note, we work with theta being the
      // angle with respect to the xy plane, which is not
      // the convention for other quadratures.  That's because
      // xi will be used in scaling flight lengths.
      d_mu[angle]     = cos(phi[a]) * d_polar->sin_theta(p);
      d_eta[angle]    = sin(phi[a]) * d_polar->sin_theta(p);
      d_xi[angle]     = d_polar->cos_theta(p);
      // Approximation for now.  It would be better to adjust
      // weights based on adjusted azimuths.
      d_weight[angle] = 2.0 * d_polar->weight(p) * delta;

    } // end polar

  } // end azimuth

  for (int a = 1; a < d_number_azimuths_octant; a++)
  {
    Insist(!soft_equiv(phi[a], phi[a-1]),
      "Repeated azimuths.  Increase space parameter or try an odd number.");
  }

  // Calculate intercepts on a square.
  d_enter.resize(2*d_number_azimuths_octant, vec_point(num_space));
  d_exit.resize(2*d_number_azimuths_octant,  vec_point(num_space));
  d_number_enter.resize(2*d_number_azimuths_octant, vec_int(2, 0));
  d_number_exit.resize(2*d_number_azimuths_octant, vec_int(2, 0));

  // First quadrant
  if (db) std::cout << std::endl;
  for (int a = 0; a < d_number_azimuths_octant; a++)
  {
    double dx = 1.0 / num_x[a];
    double dy = 1.0 / num_y[a];
    d_phi[a]  = phi[a];
    d_cos_phi[a] = cos(phi[a]);
    d_sin_phi[a] = sin(phi[a]);
    d_spacing[a] = d_sin_phi[a] * std::min(dx, 0.5);
    Assert(d_spacing[a] > 0.0);

    for (int i = 0; i < num_y[a]; i++)
    {
      d_enter[a][i]            = Point(0.0, 1.0 - (i + 0.5)*dy);
      d_exit[a][i + num_x[a]]  = Point(1.0, 1.0 - (i + 0.5)*dy);
    }
    for (int i = 0; i < num_x[a]; i++)
    {
      d_enter[a][i + num_y[a]] = Point((0.5 + i)*dx, 0.0);
      d_exit[a][i]             = Point((0.5 + i)*dx, 1.0);
    }
    d_number_enter[a][0] = num_x[a];
    d_number_enter[a][1] = num_y[a];
    d_number_exit[a][0]  = num_y[a];
    d_number_exit[a][1]  = num_x[a];
  }
  // Second quadrant, placed as mirror image.
  for (int a = 0; a < d_number_azimuths_octant;  a++)
  {
    int a2 = a + d_number_azimuths_octant;

    double dx = 1.0 / num_x[a];
    double dy = 1.0 / num_y[a];
    d_phi[a2]     = pi - phi[a];
    d_cos_phi[a2] = cos(d_phi[a2]);
    d_sin_phi[a2] = sin(d_phi[a2]);
    d_spacing[a2] = sin(phi[a]) * std::min(dx, 0.5);
    Assert(d_spacing[a2] > 0.0);

    for (int i = 0; i < num_y[a]; i++)
    {
      d_enter[a2][i]            = Point(1.0, 1.0 - (i + 0.5)*dy);
      d_exit[a2][i + num_x[a]]  = Point(0.0, 1.0 - (i + 0.5)*dy);
    }
    for (int i = 0; i < num_x[a]; i++)
    {
      d_enter[a2][i + num_y[a]] = Point(1.0 - (0.5 + i)*dx, 0.0);
      d_exit[a2][i]             = Point(1.0 - (0.5 + i)*dx, 1.0);
    }
    d_number_enter[a2][0] = num_x[a];
    d_number_enter[a2][1] = num_y[a];
    d_number_exit[a2][0]  = num_y[a];
    d_number_exit[a2][1]  = num_x[a];
    d_azimuth_weight[a2]  = d_azimuth_weight[a];
  }

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Uniform.cc
//---------------------------------------------------------------------------//
