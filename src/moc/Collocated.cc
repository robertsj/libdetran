//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Collocated.cc
 * \brief  Collocated 
 * \author Jeremy Roberts
 * \date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "Collocated.hh"

// System
#include <cmath>
#include <iostream>

namespace detran
{

Collocated::Collocated(int dim,
                       int num_azimuths_octant,
                       int multiplier,
                       int num_polar,
                       std::string polar)
  : QuadratureMOC(dim,
                  num_azimuths_octant,
                  num_polar,
                  "COLLOCATED",
                  polar)
{

  using std::cout;
  using std::endl;

  bool db = false;

  Insist(   (d_number_azimuths_octant == 3)
         or (d_number_azimuths_octant == 5)
         or (d_number_azimuths_octant == 7)
         or (d_number_azimuths_octant == 9),
         "Collocated quadrature implemented for 3, 5, 7, or 9 azimuths.");

  Insist(multiplier > 0, "The multiplier must be positive");

  int num_space = std::pow(3.0, double(d_number_azimuths_octant / 2));
  Ensure(   num_space ==  3
         or num_space ==  9
         or num_space == 27
         or num_space == 81);
  // num_space is the maximum number of tracks from any side

  if (db)
  {
    cout << "NUMBER AZIMUTHS = " << d_number_azimuths_octant << endl;
    cout << "   NUMBER SPACE = " << num_space << endl;
  }

  vec_int num_x(d_number_azimuths_octant, num_space);
  vec_int num_y(d_number_azimuths_octant, num_space);
  vec_dbl phi(d_number_azimuths_octant, 0.0);
  vec_int num_tracks(d_number_azimuths_octant, 0);

  for (int i = 0; i < d_number_azimuths_octant / 2; i++)
  {
    int a = d_number_azimuths_octant - i - 1;
    num_y[a] = std::pow(3, i); // 3^(i-1);
    num_x[i] = std::pow(3, i);
  }

  //-------------------------------------------------------------------------//
  // DEFINE QUADRATURE POINTS
  //-------------------------------------------------------------------------//

  // All quadratures are subsets of these points.
  double points[] = {0.012345051844224,
                     0.037020115873930,
                     0.110657221173896,
                     0.321750554396642,
                     0.785398163397448,
                     1.249045772398254,
                     1.460139105621001,
                     1.533776210920966,
                     1.558451274950672};

  int start;
  if (d_number_azimuths_octant == 3)
  {
    start = 3;
  }
  else if (d_number_azimuths_octant == 5)
  {
    start = 2;
  }
  else if (d_number_azimuths_octant == 7)
  {
    start = 1;
  }
  else
  {
    start = 0;
  }

  for (int a = 0; a < d_number_azimuths_octant; a++)
  {
    d_phi[a] = points[start + a];
    d_phi[2*d_number_azimuths_octant - a - 1] = pi - d_phi[a];
  }

  //-------------------------------------------------------------------------//
  // DEFINE QUADRATURE WEIGHTS
  //-------------------------------------------------------------------------//

  bool weight_flag = true;

  if (weight_flag) // Use the interpolating weights
  {
    int a = 0;
    if (d_number_azimuths_octant == 3)
    {
      double w[] = {0.751228992310453, 0.068338342173990};
      for (a = 0; a < d_number_azimuths_octant/2; a++)
      {
        cout << " weight = " << w[a] << endl;
        d_azimuth_weight[a] = w[a];
      }
      cout << " aa == " << a << endl;
      d_azimuth_weight[a] = w[a];
    }
    else if (d_number_azimuths_octant == 5)
    {
      double w[] = {0.229000262143729, 0.266237570667363, 0.580320661172713};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
    else if (d_number_azimuths_octant == 7)
    {
      double w[] = {0.102685323921243, 0.047592024817433, 0.382904372429055,
                    0.504432884459434};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
    else
    {
      double w[] = {-0.271114982393631, 0.538350840944151, -0.155409051121381,
                    0.431467880336857, 0.484206951262897};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
  }
  else // Arc length weight
  {
    int a = 0;
    if (d_number_azimuths_octant == 3)
    {
      double w[] = {0.553574358897045, 0.463647609000806};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
    else if (d_number_azimuths_octant == 5)
    {
      double w[] = {0.216203887785269, 0.337370471111776, 0.463647609000806};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
    else if (d_number_azimuths_octant == 7)
    {
      double w[] = {0.073838668523913, 0.142365219261356, 0.337370471111776,
                    0.463647609000806};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
    else
    {
      double w[] = {0.024682583859077, 0.049156084664836, 0.142365219261356,
                    0.337370471111776, 0.463647609000806};
      for (a = 0; a < d_number_azimuths_octant - 1; a++)
        d_azimuth_weight[a] = w[a];
      d_azimuth_weight[a] = w[a];
    }
  }
  for (int a = 0; a < d_number_azimuths_octant/2; a++)
  {
    int ra = d_number_azimuths_octant - a - 1;
    d_azimuth_weight[ra] =  d_azimuth_weight[a];
  }
  for (int a = 0; a < d_number_azimuths_octant; a++)
  {
    d_azimuth_weight[2*d_number_azimuths_octant - a - 1] = d_azimuth_weight[a];
  }

  //-------------------------------------------------------------------------//
  // FINALIZE QUADRATURE
  //-------------------------------------------------------------------------/

  for (int a = 0; a < d_number_azimuths_octant; a++)
  {

    for (int p = 0; p < d_number_polar; p++)
    {
      // Cardinal Index
      int angle = a * d_number_polar + p;
      d_mu[angle]     = cos(phi[a]) * d_polar->sin_theta(p);
      d_eta[angle]    = sin(phi[a]) * d_polar->sin_theta(p);
      d_xi[angle]     = d_polar->cos_theta(p);
      d_weight[angle] = 2.0 * d_polar->weight(p) * d_azimuth_weight[a];

    } // end polar
  } // azimuth

  if (db)
  {
    for (int a = 0; a < 2*d_number_azimuths_octant; a++)
    {
      int aa = a;
      if (a >= d_number_azimuths_octant)
      {
        aa = 2*d_number_azimuths_octant - a - 1;
      }
      cout << " aa = " << aa << endl;
      cout << "AZIMUTH: " <<  a << endl;
      cout << "       phi = " << d_phi[a] << endl;
      cout << "    weight = " << d_azimuth_weight[a] << endl;
      cout << "  number x = " << num_x[aa] << endl;
      cout << "  number y = " << num_y[aa] << endl;
      cout << "  number   = " << num_x[aa] + num_y[aa] << endl;
    }
  }

  //-------------------------------------------------------------------------//
  // SETUP TRACK POINTS
  //-------------------------------------------------------------------------//

  // Calculate intercepts on a square.
  d_enter.resize(2*d_number_azimuths_octant);
  d_exit.resize(2*d_number_azimuths_octant);
  d_number_enter.resize(2*d_number_azimuths_octant, vec_int(2, 0));
  d_number_exit.resize(2*d_number_azimuths_octant, vec_int(2, 0));

  // First quadrant
  for (int a = 0; a < d_number_azimuths_octant; a++)
  {
    // Number of tracks
    int nx = num_x[a];
    int ny = num_y[a];
    int n  = nx + ny;

    // Resize end points.
    d_enter[a].resize(n, 0);
    d_exit[a].resize(n, 0);

    // Horizontal and vertical steps
    double dx = 1.0 / nx;
    double dy = 1.0 / ny;

    // Fill cos and sin
    d_cos_phi[a] = cos(d_phi[a]);
    d_sin_phi[a] = sin(d_phi[a]);

    // Perpendicular distance between tracks
    d_spacing[a] = d_sin_phi[a] * dx;
    Assert(d_spacing[a] > 0.0);

    for (int i = 0; i < num_y[a]; i++)
    {
      d_enter[a][i]              = Point(0.0, 1.0 - (i + 0.5)*dy);
      d_exit[a][i + num_x[a]] = Point(1.0, 1.0 - (i + 0.5)*dy);
    }
    for (int i = 0; i < num_x[a]; i++)
    {
      d_enter[a][i + num_y[a]] = Point((0.5 + i)*dx, 0.0);
      d_exit[a][i]                = Point((0.5 + i)*dx, 1.0);
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

    // Number of tracks
    int nx = num_x[a];
    int ny = num_y[a];
    int n  = nx + ny;

    // Resize end points.
    d_enter[a2].resize(n, 0);
    d_exit[a2].resize(n, 0);

    // Horizontal and vertical steps
    double dx = 1.0 / nx;
    double dy = 1.0 / ny;

    d_cos_phi[a2] = cos(d_phi[a2]);
    d_sin_phi[a2] = sin(d_phi[a2]);
    d_spacing[a2] = d_spacing[a];
    Assert(d_spacing[a2] > 0.0);

    for (int i = 0; i < num_y[a]; i++)
    {
      d_enter[a2][i]           = Point(1.0, 1.0 - (i + 0.5)*dy);
      d_exit[a2][i + num_x[a]] = Point(0.0, 1.0 - (i + 0.5)*dy);
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
//              end of file Collocated.cc
//---------------------------------------------------------------------------//
