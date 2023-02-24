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
{
  Require(d_quadrature);
}


//----------------------------------------------------------------------------//
void TrackDB::normalize(const std::vector<double> &volume)
{
  std::vector<double> volume_appx(volume.size(), 0.0);
  //size_t na = d_quadrature->number_azimuths_octant();

  // compute approximate volumes from track widths, segment lengths, and
  // angular weights.
  for (auto &track: d_tracks)
  {

  }



//  for (size_t a = 0; a < d_tracks.size(); ++a)
//  {
//    double a_wt = d_quadrature->azimuth_weight(a % na) / detran_utilities::pi;
//    for (size_t p = 0; p < d_tracks[a].size(); ++p)
//    {
//      if (d_dimension == 3) a_wt *= d_quadrature->polar_weight(p)/2.0;
//      for (size_t t = 0; t < d_tracks[a][p].size(); t++)
//      {
//        SP_track trk =  d_tracks[a][p][t];
//        for (size_t s = 0; s < trk->number_segments(); ++s)
//        {
//          size_t region = trk->segment(s).region();
//          Assert(region < volume.size());
//          volume_appx[region] += trk->segment(s).length() * trk->width() * a_wt;
//        }
//      }
//    }
//  }
////  for (int r = 0; r < d_number_regions; r++)
////  {
////    std::cout << " REGION: " << r
////              << " volume: " << volume[r]
////              << " appx: " << appx_volume[r] << std::endl;
////  }
//  for (size_t a = 0; a < d_tracks.size(); ++a)
//  {
//    for (size_t p = 0; p < d_tracks[a].size(); ++p)
//    {
//      for (size_t t = 0; t < d_tracks[a][p].size(); ++t)
//      {
//        SP_track trk =  d_tracks[a][p][t];
//        for (int s = 0; s < d_tracks[a][p][t]->number_segments(); ++s)
//        {
//          size_t r = trk->segment(s).region();
//          trk->segment(s).scale(volume[r] / volume_appx[r]);
//        }
//      }
//    }
//  }
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
//  using std::cout;
//  using std::endl;
//  cout << endl;
//  cout << "------------" << endl;
//  cout << "TrackDB data" << endl;
//  cout << "------------" << endl;
//  cout << endl << endl;
//  cout << "dimension: " << d_dimension << endl;
//
//  for (size_t a = 0; a < d_tracks.size(); ++a)
//  {
//    cout << "  azimuth = " << a << endl;
//    for (size_t p = 0; p < d_tracks[a].size(); ++p)
//    {
//      cout << "      polar = " << p << endl;
//      for (size_t t = 0; t < d_tracks[a][p].size(); ++t)
//      {
//        Assert(d_tracks[a][p][t]);
//        cout << "        track = " << t << *d_tracks[a][p][t];
//      }
//    }
//  }
}

} // end namespace detran_geometry

//---------------------------------------------------------------------------//
//              end of file TrackDB.cc
//---------------------------------------------------------------------------//
