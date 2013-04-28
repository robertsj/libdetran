//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   QuadratureMOC.cc
 * \brief  QuadratureMOC member definitions
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#include "QuadratureMOC.hh"
#include "TabuchiYamamoto.hh"

namespace detran_angle
{

QuadratureMOC::QuadratureMOC(size_t dim,
                             size_t num_azimuths_octant,
                             size_t num_polar,
                             std::string name,
                             std::string polar)
  : Quadrature(dim,
               num_polar*num_azimuths_octant*int(std::pow(2.0, dim)),
               name+"-"+polar)
  , d_number_azimuths_octant(num_azimuths_octant)
  , d_number_polar(num_polar)
  , d_phi(2*num_azimuths_octant, 0.0)
  , d_cos_phi(2*num_azimuths_octant, 0.0)
  , d_sin_phi(2*num_azimuths_octant, 0.0)
  , d_spacing(2*num_azimuths_octant, 0.0)
  , d_azimuth_weight(2*num_azimuths_octant, 0.0)
{

  // Support for 2D right now
  Require(dim == 2);

  // Construct the polar quadrature.
  if (polar == "TY")
  {
    d_polar = new TabuchiYamamoto(num_polar);
  }
  else
  {
    THROW("Unsupported polar quadrature specified:" + polar);
  }

}

void QuadratureMOC::display_tracks() const
{
  using std::cout;
  using std::endl;

  cout << "MOC Quadrature: " << d_name << endl;
  cout << endl;
  cout << "Unscaled track data: " << endl;
  for (size_t a = 0; a < 2*d_number_azimuths_octant; ++a)
  {
    cout << "*** azimuth = " << a << endl;
    cout << "        phi = " << d_phi[a] << endl;
    cout << "    cos_phi = " << d_cos_phi[a] << endl;
    cout << "    sin_phi = " << d_sin_phi[a] << endl;
    cout << "    spacing = " << d_spacing[a] << endl;
    for (size_t t = 0; t < d_exit[a].size(); ++t)
    {
      cout << "****** track = " << t;
      cout << " enter = " << d_enter[a][t]
           << " exit = "  << d_exit[a][t] << endl;

    }
  }

}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of file QuadratureMOC.cc
//---------------------------------------------------------------------------//
