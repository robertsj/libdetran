//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TrackDB.cc
 *  @brief TrackDB member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "TrackDB.hh"
#include <iostream>

namespace detran_geometry
{

//----------------------------------------------------------------------------//
TrackDB::TrackDB(SP_quadrature quad,
                 const size_t  num_regions)
  : d_quadrature(quad)
  , d_number_regions(num_regions)
{
  Require(d_quadrature);
  Require(d_number_regions > 0);
  d_number_azimuths = d_quadrature->number_azimuths_octant() * 2;


  d_tracks.resize(d_number_azimuths);

  Ensure(d_number_azimuths > 0);
}

//----------------------------------------------------------------------------//
TrackDB::SP_track TrackDB::track(const size_t a, const size_t t)
{
  Require(a < d_tracks.size());
  Require(t < d_tracks[a].size());
  return d_tracks[a][t];
}

//----------------------------------------------------------------------------//
TrackDB::size_t TrackDB::number_tracks_angle(const size_t a) const
{
  Require(a < d_tracks.size());
  return d_tracks[a].size();
}

//----------------------------------------------------------------------------//
TrackDB::size_t TrackDB::number_angles() const
{
  return d_tracks.size();
}

//----------------------------------------------------------------------------//
void TrackDB::add_track(const size_t a, SP_track t)
{
  Requirev(a < d_tracks.size(), AsString(a)+" !< "+AsString(d_tracks.size()));
  Require(t);
  d_tracks[a].push_back(t);
}

//----------------------------------------------------------------------------//
void TrackDB::setup_angle(const size_t a,
                          const double c_phi,
                          const double s_phi)
{
  Require(a < d_tracks.size());
  d_cos_phi[a] = c_phi;
  d_sin_phi[a] = s_phi;
}

//----------------------------------------------------------------------------//
void TrackDB::normalize(const vec_dbl &volume)
{
  vec_dbl appx_volume(volume.size(), 0.0);
  for (size_t a = 0; a < d_tracks.size(); a++)
  {
    for (size_t t = 0; t < d_tracks[a].size(); t++)
    {
      for (int s = 0; s < d_tracks[a][t]->number_segments(); s++)
      {
        int region = d_tracks[a][t]->segment(s).region();
        Assert(region < volume.size());
        appx_volume[region] +=
          d_tracks[a][t]->segment(s).length() * d_tracks[a][t]->width() *
            d_quadrature->azimuth_weight(a) / detran_utilities::pi;
      }
    }
  }
//  for (int r = 0; r < d_number_regions; r++)
//  {
//    std::cout << " REGION: " << r
//              << " volume: " << volume[r]
//              << " appx: " << appx_volume[r] << std::endl;
//  }
  for (size_t a = 0; a < d_tracks.size(); a++)
  {
    for (size_t t = 0; t < d_tracks[a].size(); t++)
    {
      for (int s = 0; s < d_tracks[a][t]->number_segments(); s++)
      {
        int r = d_tracks[a][t]->segment(s).region();
        d_tracks[a][t]->segment(s).scale(volume[r]/appx_volume[r]);
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
  cout << endl;
  for (size_t a = 0; a < d_tracks.size(); a++)
  {
    cout << "    azimuth = " << a << endl;
    for (size_t t = 0; t < d_tracks[a].size(); t++)
    {
      cout << "       track = " << t << *d_tracks[a][t] << endl;
    }
  }
}

} // end namespace detran_geometry

//---------------------------------------------------------------------------//
//              end of file TrackDB.cc
//---------------------------------------------------------------------------//
