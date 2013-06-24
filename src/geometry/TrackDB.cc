//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TrackDB.cc
 *  @brief TrackDB member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "TrackDB.hh"
#include "utilities/TinyVector.hh"
#include <iostream>

namespace detran_geometry
{
#define COUT(c) std::cout << c << std::endl;
//----------------------------------------------------------------------------//
TrackDB::TrackDB(SP_quadrature q)
  : d_quadrature(q)
  , d_number_polar(1)
{
  Require(d_quadrature);

  d_dimension = d_quadrature->dimension();
  if (d_dimension == 2)
  {
    d_number_azimuths = d_quadrature->number_azimuths_octant() * 2;
  }
  else
  {
    d_number_azimuths = d_quadrature->number_azimuths_octant() * 4;
    d_number_polar    = d_quadrature->number_polar_octant();
  }

  d_tracks.resize(d_number_azimuths, vec2_track(d_number_polar));

  Ensure(d_number_azimuths > 0);
}

//----------------------------------------------------------------------------//
TrackDB::SP_track TrackDB::track(c_size_t a, c_size_t p, c_size_t t)
{
  Require(a < d_tracks.size());
  Require(p < d_tracks[a].size());
  Require(t < d_tracks[a][p].size());
  return d_tracks[a][p][t];
}

//----------------------------------------------------------------------------//
TrackDB::iterator_angle TrackDB::begin(c_size_t a, c_size_t p)
{
  Require(d_dimension = 2 ? p == 0 : true);
  size_t aa = a % d_quadrature->number_azimuths_octant();
  size_t pp = p % d_quadrature->number_polar_octant();
  Ensure(aa < d_tracks.size());
  Ensure(pp < d_tracks[a].size());
  return d_tracks[aa][pp].begin();
}

//----------------------------------------------------------------------------//
TrackDB::iterator_angle TrackDB::end(c_size_t a, c_size_t p)
{
  Require(d_dimension = 2 ? p == 0 : true);
  size_t aa = a % d_quadrature->number_azimuths_octant();
  size_t pp = p % d_quadrature->number_polar_octant();
  Ensure(aa < d_tracks.size());
  Ensure(pp < d_tracks[a].size());
  return d_tracks[aa][pp].end();
}

//----------------------------------------------------------------------------//
TrackDB::size_t TrackDB::number_tracks(c_size_t a, c_size_t p) const
{
  Require(a < d_tracks.size());
  Require(p < d_tracks[a].size());
  return d_tracks[a][p].size();
}

//----------------------------------------------------------------------------//
void TrackDB::add_track(c_size_t a, c_size_t p, SP_track t)
{
  //COUT("a=" << a << " p=" << p)
  Requirev(a < d_tracks.size(), AsString(a)+" !< "+AsString(d_tracks.size()));
  Require(p < d_tracks[a].size());
  Require(t);
  d_tracks[a][p].push_back(t);
}

//----------------------------------------------------------------------------//
void TrackDB::normalize(const vec_dbl &volume)
{
  vec_dbl volume_appx(volume.size(), 0.0);
  size_t na = d_quadrature->number_azimuths_octant();
  for (size_t a = 0; a < d_tracks.size(); ++a)
  {
    double a_wt = d_quadrature->azimuth_weight(a % na) / detran_utilities::pi;
    for (size_t p = 0; p < d_tracks[a].size(); ++p)
    {
      if (d_dimension == 3) a_wt *= d_quadrature->polar_weight(p)/2.0;
      for (size_t t = 0; t < d_tracks[a][p].size(); t++)
      {
        SP_track trk =  d_tracks[a][p][t];
        for (size_t s = 0; s < trk->number_segments(); ++s)
        {
          size_t region = trk->segment(s).region();
          Assert(region < volume.size());
          volume_appx[region] += trk->segment(s).length() * trk->width() * a_wt;
        }
      }
    }
  }
//  for (int r = 0; r < d_number_regions; r++)
//  {
//    std::cout << " REGION: " << r
//              << " volume: " << volume[r]
//              << " appx: " << appx_volume[r] << std::endl;
//  }
  for (size_t a = 0; a < d_tracks.size(); ++a)
  {
    for (size_t p = 0; p < d_tracks[a].size(); ++p)
    {
      for (size_t t = 0; t < d_tracks[a][p].size(); ++t)
      {
        SP_track trk =  d_tracks[a][p][t];
        for (int s = 0; s < d_tracks[a][p][t]->number_segments(); ++s)
        {
          size_t r = trk->segment(s).region();
          trk->segment(s).scale(volume[r] / volume_appx[r]);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------//
TrackDB::size_t TrackDB::number_azimuths() const
{
  return d_number_azimuths;
}

//----------------------------------------------------------------------------//
TrackDB::size_t TrackDB::number_polar() const
{
  return d_number_polar;
}

//----------------------------------------------------------------------------//
TrackDB::size_t TrackDB::dimension() const
{
  return d_dimension;
}


//----------------------------------------------------------------------------//
struct TrackEntranceCompare
{
  typedef detran_utilities::TinyVector<int, 3> TV;
  TrackEntranceCompare(const int quadrant)
  {
    const int X=1, Y=2, Z=3;
    if      (quadrant == 0)
      d_order = TV( X,-Y, Z);
    else if (quadrant == 1)
      d_order = TV( X, Y, Z);
    else if (quadrant == 2)
      d_order = TV( Y,-X, Z);
    else if (quadrant == 3)
      d_order = TV(-X,-Y, Z);
    d_sign = TV(d_order[0] > 0, d_order[1] > 0, d_order[2] > 0);
  }
  bool operator()(const Track::SP_track a, const Track::SP_track b)
  {
    //COUT("lala")
    Point Pa = a->enter();
    Point Pb = b->enter();
    for (int d = 0; d < 3; ++d)
    {

      int  i = std::abs(d_order[d])-1;
      //COUT("tv  d = " << d << " i=" << i)
      bool L = Pa[i] < Pb[i];
      bool G = Pa[i] > Pb[i];
      if      (d_sign[d] ? L : G)
      {
        return true;
      }
      else if (d_sign[d] ? G : L)
      {
         return false;
      }
    }
    return false;
  }
  TV d_order;
  TV d_sign;
};

//----------------------------------------------------------------------------//
void TrackDB::sort()
{
  size_t nq = (d_dimension == 2) ? 2 : 4;
  size_t na = d_quadrature->number_azimuths_octant();
  for (size_t q = 0; q < nq; ++q)
  {
    for (size_t a = 0; a < na; ++a)
    {
      for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
      {
        std::sort(d_tracks[a+na*q][p].begin(),
                  d_tracks[a+na*q][p].end(),
                  TrackEntranceCompare(q));
      }
    }
  }
}

//----------------------------------------------------------------------------//
void TrackDB::display() const
{
  using std::cout;
  using std::endl;
  cout << endl;
  cout << "------------" << endl;
  cout << "TrackDB data" << endl;
  cout << "------------" << endl;
  cout << endl << endl;
  cout << "dimension: " << d_dimension << endl;

  for (size_t a = 0; a < d_tracks.size(); ++a)
  {
    cout << "  azimuth = " << a << endl;
    for (size_t p = 0; p < d_tracks[a].size(); ++p)
    {
      cout << "      polar = " << p << endl;
      for (size_t t = 0; t < d_tracks[a][p].size(); ++t)
      {
        Assert(d_tracks[a][p][t]);
        cout << "        track = " << t << *d_tracks[a][p][t];
      }
    }
  }
}

} // end namespace detran_geometry

//---------------------------------------------------------------------------//
//              end of file TrackDB.cc
//---------------------------------------------------------------------------//
