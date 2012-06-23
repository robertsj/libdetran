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
