//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   Quadrature.cc
 *  @author Jeremy Roberts
 *  @date   Mar 23, 2012
 *  @brief  Quadrature class member definitions
 */
//----------------------------------------------------------------------------//

#include "angle/Quadrature.hh"
#include <iostream>
#include <cstdio>
#include <cmath>

namespace detran_angle
{

//----------------------------------------------------------------------------//
Quadrature::Quadrature(const size_t      dim,
                       const size_t      number_angles,
                       const std::string name)
  : d_dimension(dim)
  , d_number_angles(number_angles)
  , d_number_octants(std::pow((float)2, (int)dim))
  , d_number_angles_octant(number_angles/d_number_octants)
  , d_weight(d_number_angles_octant, 0.0)
  , d_mu(d_number_angles_octant, 0.0)
  , d_eta(d_number_angles_octant, 0.0)
  , d_xi(d_number_angles_octant, 0.0)
  , d_name(name)
  , d_octant_sign(8, vec_dbl(3, 0.0))
  , d_incident_octants(2*dim, vec_int(d_number_octants/2))
  , d_outgoing_octants(2*dim, vec_int(d_number_octants/2))
  , d_adjoint(false)
{
  // Preconditions
  Insist((dim > 0) && (dim < 4),
    "The quadrature dimension MUST be 1, 2, or 3.");

  // Define the signs for all eight octants.
  // first
  d_octant_sign[0][0] =  1.0;
  d_octant_sign[0][1] =  1.0;
  d_octant_sign[0][2] =  1.0;
  // second
  d_octant_sign[1][0] = -1.0;
  d_octant_sign[1][1] =  1.0;
  d_octant_sign[1][2] =  1.0;
  // third
  d_octant_sign[2][0] = -1.0;
  d_octant_sign[2][1] = -1.0;
  d_octant_sign[2][2] =  1.0;
  // fourth
  d_octant_sign[3][0] =  1.0;
  d_octant_sign[3][1] = -1.0;
  d_octant_sign[3][2] =  1.0;
  // fifth
  d_octant_sign[4][0] =  1.0;
  d_octant_sign[4][1] =  1.0;
  d_octant_sign[4][2] = -1.0;
  // sixth
  d_octant_sign[5][0] = -1.0;
  d_octant_sign[5][1] =  1.0;
  d_octant_sign[5][2] = -1.0;
  // seventh
  d_octant_sign[6][0] = -1.0;
  d_octant_sign[6][1] = -1.0;
  d_octant_sign[6][2] = -1.0;
  // eighth
  d_octant_sign[7][0] =  1.0;
  d_octant_sign[7][1] = -1.0;
  d_octant_sign[7][2] = -1.0;

  // Set incident and outgoing
  if (d_dimension == 1)
  {
    d_incident_octants[0][0] = 0;
    d_incident_octants[1][0] = 1;
    d_outgoing_octants[0][0] = 1;
    d_outgoing_octants[1][0] = 0;
  }
  else if (d_dimension == 2)
  {
    int inc[4][2] = {{3, 0}, {1, 2},
                     {0, 1}, {2, 3}};
    int out[4][2] = {{1, 2}, {3, 0},
                     {2, 3}, {0, 1}};
//    int inc[4][2] = {{0, 3}, {2, 1},
//                     {1, 0}, {3, 2}};
//    int out[4][2] = {{2, 1}, {0, 3},
//                     {3, 2}, {1, 0}};
    for (int s = 0; s < 4; ++s)
    {
      for (int o = 0; o < 2; ++o)
      {
        d_incident_octants[s][o] = inc[s][o];
        d_outgoing_octants[s][o] = out[s][o];
      }
    }
  }
  else if (d_dimension == 3)
  {
    // vertical surface: r->l, b->t
    // horizontal      : counter clockwise
    int inc[6][4] = {{7,4,3,0}, {5,6,1,2},
                     {4,5,0,1}, {6,7,2,3},
                     {0,1,2,3}, {4,5,6,7}};
    int out[6][4] = {{5,6,1,2}, {7,4,3,0},
                     {6,7,2,3}, {4,5,0,1},
                     {4,5,6,7}, {0,1,2,3}};
    for (int s = 0; s < 6; ++s)
    {
      for (int o = 0; o < 4; ++o)
      {
        d_incident_octants[s][o] = inc[s][o];
        d_outgoing_octants[s][o] = out[s][o];
      }
    }
  }

}

//----------------------------------------------------------------------------//
Quadrature::~Quadrature()
{
  /* ... */
}

//----------------------------------------------------------------------------//
void Quadrature::set_adjoint(const bool v)
{

  // If our adjoint flag is unchanged, return. Otherwise reset flag
  // and reset the octant multipliers by multiplying by -1.
  if (v == d_adjoint)
    return;
  else
    d_adjoint = v;

  for (int o = 0; o < 8; o++)
    for (int d = 0; d < 3; d++)
      d_octant_sign[o][d] *= -1.0;

}

//----------------------------------------------------------------------------//
const Quadrature::vec_int&
Quadrature::incident_octant(const size_t s)
{
  Require(s < d_incident_octants.size());
  return d_incident_octants[s];
}

//----------------------------------------------------------------------------//
const Quadrature::vec_int&
Quadrature::outgoing_octant(const size_t s)
{
  Require(s < d_outgoing_octants.size());
  return d_outgoing_octants[s];
}

//----------------------------------------------------------------------------//
void Quadrature::display() const
{

  using std::cout;
  using std::endl;
  using std::printf;

  cout << endl;
  cout << d_name << " abscissa and weights: " << endl << endl;

  double weight_sum = 0;

  if (d_dimension == 1)
  {
    cout << "   m            mu                  wt       " << endl;
    cout << "  ---   ------------------  -----------------" << endl;
    for (size_t ix = 0; ix < d_number_angles_octant; ++ix)
    {
      printf ("%4i    %16.13f   %16.13f   \n", ix, d_mu[ix], d_weight[ix] );
      weight_sum += d_weight[ix];
    }
  }
  else
  {
    cout << "   m            mu                 eta                xi                 wt       " << endl;
    cout << "  ---   -----------------  -----------------  -----------------  -----------------" << endl;
    for (size_t ix = 0; ix < d_number_angles_octant; ++ix)
    {
      printf ("%4i    %16.13f   %16.13f   %16.13f   %16.13f   \n",
              ix, d_mu[ix], d_eta[ix], d_xi[ix], d_weight[ix] );
      weight_sum += d_weight[ix];
    }
  }

  cout << endl << "  The sum of the weights is " << weight_sum << endl << endl;
}

} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of Quadrature.cc
//----------------------------------------------------------------------------//
