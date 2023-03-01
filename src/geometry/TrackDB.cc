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
TrackDB::TrackDB(SP_quadrature q, const map_track& tracks)
  : d_quadrature(q)
  , d_tracks(tracks)
{
  Require(d_quadrature);
}

//----------------------------------------------------------------------------//
std::vector<double> TrackDB::volume() const
{
  // compute volumes from track widths, segment lengths, and
  // angular weights.
  std::map<int, double> volume_appx;
  double total_angle_weight = 0;
  for (auto &angle_tracks: d_tracks)
  {
    int a = angle_tracks.first.second;
    double angle_weight = d_quadrature->weight(a);
    total_angle_weight += angle_weight;
    for (auto track: angle_tracks.second)
    {
      for (auto segment: track->segments())
      {
        int region = segment.region();
        double dV = segment.length() * track->spatial_weight() * angle_weight;
        COUT(" region = " << region << " dV = " << dV << " angle_weight = " << angle_weight);
        if (volume_appx.count(region)==0)
          volume_appx[region] = 0.0;
        volume_appx[region] += dV;
      }
    }
  }
  COUT(" total_angle_weight = " << total_angle_weight);
  // could check angular weight
  std::vector<double> volume_appx_v(volume_appx.size(), 0.0);
  for (auto &v: volume_appx)
  {
    auto r = v.first;
    Assertv(r < volume_appx.size(), "missing region?"); // regions should be 0, 1, ...
    volume_appx_v[r] = v.second / total_angle_weight;
  }

  return volume_appx_v;
}


//----------------------------------------------------------------------------//
void TrackDB::normalize(const std::vector<double> &ref_volume)
{
  std::vector<double> appx_volume = volume();
  Assert(ref_volume.size()==appx_volume.size());

  for (int r = 0; r < ref_volume.size(); r++)
  {
    std::cout << " REGION: " << r
              << " volume: " << ref_volume[r]
              << "   appx: " << appx_volume[r] << std::endl;
  }

  for (auto &angle_tracks: d_tracks)
   {
     for (auto track: angle_tracks.second)
     {
       for (auto& segment: track->segments())
       {
         int r = segment.region();
         segment.scale(ref_volume[r] / appx_volume[r]);
       }
     }
   }
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
  bool operator()(Track::SP_track a, Track::SP_track b)
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
  // \todo Ensure this sort is needed or delete.
//  size_t nq = (d_dimension == 2) ? 2 : 4;
//  size_t na = d_quadrature->number_azimuths_octant();
//  for (size_t q = 0; q < nq; ++q)
//  {
//    for (size_t a = 0; a < na; ++a)
//    {
//      for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
//      {
//        std::sort(d_tracks[a+na*q][p].begin(),
//                  d_tracks[a+na*q][p].end(),
//                  TrackEntranceCompare(q));
//      }
//    }
//  }
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
  for (auto &angle_tracks : d_tracks)
  {
    int o = angle_tracks.first.first;
    int a = angle_tracks.first.second;
    cout << "  octant =  " << o << "  angle = " << a << endl;
    for (auto track : angle_tracks.second)
    {
      cout << *track << endl;
    }
  }
}

} // end namespace detran_geometry

//---------------------------------------------------------------------------//
//              end of file TrackDB.cc
//---------------------------------------------------------------------------//
