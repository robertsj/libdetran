//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.cc
 *  @brief ProductQuadrature
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "angle/ProductQuadrature.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
ProductQuadrature::ProductQuadrature(const size_t       dim,
                                     const size_t       na,
                                     const size_t       np,
                                     const std::string &name,
                                     const bool         normalize)
  : Quadrature(dim,
               na * np * std::pow((float)2, (int)dim),
               name)
  , d_number_azimuth_octant(na)
  , d_number_polar_octant(np)
  , d_phi(d_number_azimuth_octant, 0.0)
  , d_cos_phi(d_number_azimuth_octant, 0.0)
  , d_sin_phi(d_number_azimuth_octant, 0.0)
  , d_azimuth_weight(d_number_azimuth_octant, 0.0)
  , d_cos_theta(d_number_polar_octant, 0.0)
  , d_sin_theta(d_number_polar_octant, 0.0)
  , d_polar_weight(d_number_polar_octant, 0.0)
  , d_number_polar(dim * 2)
  , d_number_azimuths(dim * 2)
  , d_incident_indices(2 * dim)
  , d_outgoing_indices(2 * dim)
{
  Insist(dim > 1, "Product quadratures make no sense in one dimension");

  for (size_t s = 0; s < 2 * dim; ++s)
  {

    if (s < 4)
    {
      // vertical
      d_number_polar[s]   = dim == 2 ? np : 2 * np;
      d_number_azimuths[s] = 2 * na;
    }
    else
    {
      // horizontal
      d_number_polar[s]   = np;
      d_number_azimuths[s] = 4 * na;
    }
    d_incident_indices[s].resize(d_number_azimuths[s],
                                 vec1_index_t(d_number_polar[s]));
    d_outgoing_indices[s].resize(d_number_azimuths[s],
                                 vec1_index_t(d_number_polar[s]));
    for (size_t oo = 0; oo < (dim - 1) * 2; ++oo)
    {
      size_t o_inc = incident_octant(s)[oo];
      size_t o_out = outgoing_octant(s)[oo];
      for (size_t a = 0; a < number_angles_octant(); ++a)
      {
        size_t pol_inc = this->polar(a);
        size_t pol_out = this->polar(a);
        size_t azi_inc = this->azimuth(a);
        size_t azi_out = this->azimuth(a);
        if (s < 4)
        {
          // vertical
          if (dim == 3)
          {
            pol_inc = o_inc <= 3 ? number_polar_octant() - pol_inc - 1
                                 : pol_inc + number_polar_octant();
            pol_out = o_out <= 3 ? number_polar_octant() - pol_out - 1
                                 : pol_out + number_polar_octant();
          }
          // Right-to-left, 0 to pi for phi.
          if (s == 0)
          {
            // WEST
            if (o_inc == 3 || o_inc == 7)
              azi_inc = number_azimuths_octant() - azi_inc - 1;
            else
              azi_inc += number_azimuths_octant();
            if (o_out == 1 || o_out == 5)
              azi_out = number_azimuths_octant() - azi_out - 1;
            else
              azi_out += number_azimuths_octant();
          }
          else if (s == 1)
          {
            // EAST
            if (o_inc == 1 || o_inc == 5)
              azi_inc = number_azimuths_octant() - azi_inc - 1;
            else
              azi_inc += number_azimuths_octant();
            if (o_out == 3 || o_out == 7)
              azi_out = number_azimuths_octant() - azi_out - 1;
            else
              azi_out += number_azimuths_octant();
          }
          else if (s == 2)
          {
            // SOUTH
            if (o_inc == 1 || o_inc == 5)
              azi_inc = 2 * number_azimuths_octant() - azi_inc - 1;
            else
              azi_inc = azi_inc * 1;
            if (o_out == 3 || o_out == 7)
              azi_out = 2 * number_azimuths_octant() - azi_out - 1;
            else
              azi_out = azi_inc * 1;
          }
          else
          {
            // NORTH
            if (o_inc == 3 || o_inc == 7)
              azi_inc = 2 * number_azimuths_octant() - azi_inc - 1;
            else
              azi_inc = azi_inc * 1;
            if (o_out == 1 || o_out == 5)
              azi_out = 2 * number_azimuths_octant() - azi_out - 1;
            else
              azi_out = azi_out * 1;
          }
        }
        else
        {
          // horizontal
          if (s == 4)
          {
            Assert(o_out > 3);
            Assert(o_inc < 4);
            pol_out = 1 * pol_inc;
            pol_inc = number_polar_octant() - pol_inc - 1;
          }
          else
          {
            Assert(o_inc > 3);
            Assert(o_out < 4);
            pol_out = number_polar_octant() - pol_inc - 1;
            pol_inc = 1 * pol_inc;
          }
          //
          if      (o_inc == 0 || o_inc == 4)
            azi_inc = 1 * azi_inc;
          else if (o_inc == 1 || o_inc == 5)
            azi_inc = 2 * number_azimuths_octant() - azi_inc - 1;
          else if (o_inc == 2 || o_inc == 6)
            azi_inc = 2 * number_azimuths_octant() + azi_inc;
          else if (o_inc == 3 || o_inc == 7)
            azi_inc = 4 * number_azimuths_octant() - azi_inc - 1;

          //
          if      (o_out == 0 || o_out == 4)
            azi_out = 1 * azi_out;
          else if (o_out == 1 || o_out == 5)
            azi_out = 2 * number_azimuths_octant() - azi_out - 1;
          else if (o_out == 2 || o_out == 6)
            azi_out = 2 * number_azimuths_octant() + azi_out;
          else if (o_out == 3 || o_out == 7)
            azi_out = 4 * number_azimuths_octant() - azi_out - 1;
          //printf("-- %2i %2i %2i %2i %2i %2i \n", s, oo, o_inc, o_out, azi_inc, azi_out);

        }
        d_incident_indices[s][azi_inc][pol_inc].octant = o_inc;
        d_incident_indices[s][azi_inc][pol_inc].angle  = a;
        d_incident_indices[s][azi_inc][pol_inc].io_octant = oo;
        //
        d_outgoing_indices[s][azi_out][pol_out].octant = o_out;
        d_outgoing_indices[s][azi_out][pol_out].angle  = a;
        d_outgoing_indices[s][azi_out][pol_out].io_octant = oo;
      } // angle
    } // octant
  } // side

}

void ProductQuadrature::display_indices()
{
  printf("INCIDENT INDICES\n");
  for (int s = 0; s < 2 * d_dimension; ++s)
  {
    for (int az = 0; az < this->number_azimuths(s); ++az)
    {
      for (int p = 0; p < this->number_polar(s); ++p)
      {
        int o = incident_index(s, az, p).octant;
        int a = incident_index(s, az, p).angle;
        printf("-- %2i %2i %2i %2i %2i  %8.4f  %8.4f  %8.4f \n",
                   s, o, a, az, p,
                   mu(o, a), eta(o, a), xi(o, a));
      }
    }
  }
  printf("OUTGOING INDICES\n");
  for (int s = 0; s < 2 * d_dimension; ++s)
  {
    for (int az = 0; az < this->number_azimuths(s); ++az)
    {
      for (int p = 0; p < this->number_polar(s); ++p)
      {
        int o = outgoing_index(s, az, p).octant;
        int a = outgoing_index(s, az, p).angle;
        printf("-- %2i %2i %2i %2i %2i  %8.4f  %8.4f  %8.4f \n",
                   s, o, a, az, p,
                   mu(o, a), eta(o, a), xi(o, a));
      }
    }
  }
}

//----------------------------------------------------------------------------//
ProductQuadrature::~ProductQuadrature()
{
  /* ... */
}

//----------------------------------------------------------------------------//
double ProductQuadrature::phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_phi[a];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::cos_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_cos_phi[a];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::sin_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_sin_phi[a];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::azimuth_weight(const size_t a) const
{
  Require(a < d_number_azimuth_octant);
  return d_azimuth_weight[a];
}

//----------------------------------------------------------------------------//
double ProductQuadrature::cos_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_cos_theta[p];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::sin_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_sin_theta[p];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::polar_weight(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_polar_weight[p];
}

//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::angle(const size_t a, const size_t p) const
{
  Require(a < d_number_azimuth_octant);
  Require(p < d_number_polar_octant);
  return p + a * d_number_polar_octant;
}
//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::azimuth(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = (angle % d_number_angles_octant) / d_number_polar_octant;
  Ensure(tmp < d_number_azimuth_octant);
  return tmp;
}
//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::polar(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = angle % d_number_polar_octant;
  Ensure(tmp < d_number_polar_octant);
  return tmp;
}

//----------------------------------------------------------------------------//
void ProductQuadrature::build()
{
  // Construct the product set
  double scale = 1.0;
  if (d_dimension == 2) scale = 2.0;
  double wt_tot = 0.0;
  size_t n = 0;
  double pwt = 0.0;
  for (size_t a = 0; a < number_azimuths_octant(); ++a)
  {
    pwt = 0;
    for (size_t p = 0; p < number_polar_octant(); ++p, ++n)
    {
      d_mu[n]     = d_sin_theta[p] * d_cos_phi[a];
      d_eta[n]    = d_sin_theta[p] * d_sin_phi[a];
      d_xi[n]     = d_cos_theta[p];
      d_weight[n] = scale * d_polar_weight[p] * d_azimuth_weight[a];
      wt_tot += d_weight[n]; pwt += d_polar_weight[p];
      //printf("%24.18e \n", d_sin_theta[p]);
    } // end polar loop
  } // end azimuth loop
  wt_tot *= d_number_octants;
  //std::cout << "pwt=" <<  1.0-pwt << " " << pwt << std::endl;
  verify();
//  if (!detran_utilities::soft_equiv(wt_tot, 1.0/angular_norm(d_dimension)))
//  {
//    std::cout << "***NOTICE*** sum of weights - 1/angular_norm = "
//              << wt_tot-1.0/angular_norm(d_dimension) << std::endl;
//  }
//  Ensurev(detran_utilities::soft_equiv(wt_tot, 1.0/angular_norm(d_dimension), 1e-11),
//    AsString(wt_tot) + " and not " + AsString(wt_tot-1.0/angular_norm(d_dimension)));
}

//----------------------------------------------------------------------------//
void ProductQuadrature::verify() const
{
  // verify polar ordering
  double last_p = 0.0;
  for (size_t p = 0; p < d_number_polar_octant; ++p)
  {
    Insist(d_cos_theta[p] > last_p, "Non-monotonic increasing polar cosine");
    last_p = d_cos_theta[p];
    // verify azimuth
    double last_mu = 1.0;
    for (size_t a = 0; a < d_number_azimuth_octant; ++a)
    {
      Insist(d_mu[angle(a, p)] < last_mu, "Non-monotonic increasing mu");
      last_mu = d_mu[angle(a, p)];
    }
  }
}

//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::number_azimuths(const size_t s) const
{
  Require(s < d_number_azimuths.size());
  return d_number_azimuths[s];
}

//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::number_polar(const size_t s) const
{
  Require(s < d_number_polar.size());
  return d_number_polar[s];
}
//----------------------------------------------------------------------------//
ProductQuadrature::index_t
ProductQuadrature::incident_index(const size_t s,
                                  const size_t a,
                                  const size_t p) const
{
  Require(d_number_polar.size());
  Require(a < d_number_azimuths[s]);
  Require(p < d_number_polar[s]);
  return d_incident_indices[s][a][p];
}

//----------------------------------------------------------------------------//
ProductQuadrature::index_t
ProductQuadrature::outgoing_index(const size_t s,
                                  const size_t a,
                                  const size_t p) const
{
  Require(d_number_polar.size());
  Require(a < d_number_azimuths[s]);
  Require(p < d_number_polar[s]);
  return d_outgoing_indices[s][a][p];
}

} // namespace detran_angle

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.cc
//----------------------------------------------------------------------------//
