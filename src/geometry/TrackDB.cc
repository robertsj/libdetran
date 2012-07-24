//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TrackDB.cc
 * \brief  TrackDB 
 * \author Jeremy Roberts
 * \date   Jun 23, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "TrackDB.hh"

// System
#include <iostream>

namespace detran
{

void TrackDB::normalize(vec_dbl &volume)
{
  vec_dbl appx_volume(volume.size(), 0.0);
  for (int a = 0; a < d_tracks.size(); a++)
  {
    for (int t = 0; t < d_tracks[a].size(); t++)
    {
      for (int s = 0; s < d_tracks[a][t]->number_segments(); s++)
      {
        int region = d_tracks[a][t]->segment(s).region();
        Assert(region < d_number_regions);
        appx_volume[region] +=
          d_tracks[a][t]->segment(s).length() * d_spacing[a] *
            d_quadrature->azimuth_weight(a) / pi;
      }
    }
  }
//  for (int r = 0; r < d_number_regions; r++)
//  {
//    std::cout << " REGION: " << r
//              << " volume: " << volume[r]
//              << " appx: " << appx_volume[r] << std::endl;
//  }
  for (int a = 0; a < d_tracks.size(); a++)
  {
    for (int t = 0; t < d_tracks[a].size(); t++)
    {
      for (int s = 0; s < d_tracks[a][t]->number_segments(); s++)
      {
        int r = d_tracks[a][t]->segment(s).region();
        d_tracks[a][t]->segment(s).scale(volume[r]/appx_volume[r]);
      }
    }
  }

}

void TrackDB::display() const
{
  using std::cout;
  using std::endl;
  cout << endl;
  cout << "------------" << endl;
  cout << "TrackDB data" << endl;
  cout << "------------" << endl;
  cout << endl;
  for (int a = 0; a < d_tracks.size(); a++)
  {
    cout << "    azimuth = " << a << endl;
    for (int t = 0; t < d_tracks[a].size(); t++)
    {
      cout << "       track = " << t << *d_tracks[a][t] << endl;
    }
  }

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file TrackDB.cc
//---------------------------------------------------------------------------//
