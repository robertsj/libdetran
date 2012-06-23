//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   QuadratureMOC.hh
 * \brief  QuadratureMOC class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef QUADRATUREMOC_HH_
#define QUADRATUREMOC_HH_

// Detran
#include "Quadrature.hh"
#include "PolarQuadrature.hh"
#include "Point.hh"

// Utilities
#include "Definitions.hh"

// System
#include <cmath>

namespace detran
{

/*!
 *  \class QuadratureMOC
 *  \brief Quadrature specification for MOC transport
 *
 *  In general, quadrature for MOC has intrinsic space-angle
 *  coupling.  Here, spatial aspect is incorporated
 *  by defining the start and end points of each track on
 *  a unit cell (i.e. 1x1 cell).  These points can be
 *  scaled by the tracker.
 *
 *  Also, several quantities of anticipated use in response
 *  generation are recorded.
 *
 */
class QuadratureMOC: public Quadrature
{

public:

  /// \name Useful Typedefs
  /// \{

  typedef SP<QuadratureMOC>       SP_quadrature;
  typedef Quadrature              Base;
  typedef Base::SP_quadrature     SP_base;
  typedef SP<PolarQuadrature>     SP_polar;
  typedef std::vector<Point>      vec_point;
  typedef std::vector<vec_point>  vec2_point;

  /// \}

  /*!
   *  \brief Constructor
   *  \param order                Quadrature order (\todo is this necessary?)
   *  \param dim                  Problem dimension (2 or 3d)
   *  \param num_azimuths_octant  Number azimuths per octant
   *  \param num_polar            Number of polar angles
   */
  QuadratureMOC(int dim,
                int num_azimuths_octant,
                int num_polar,
                std::string name,
                std::string polar);

  int number_azimuths_octant() const
  {
    return d_number_azimuths_octant;
  }

  int number_polar() const
  {
    return d_number_polar;
  }

  int number_tracks(u_int a)
  {
    Require(a < 2*d_number_azimuths_octant);
    return d_exit[a].size();
  }

  int number_enter(u_int a, u_int xy) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(xy == 0 or xy == 1);
    return d_number_enter[a][xy];
  }

  int number_exit(u_int a, u_int xy) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(xy == 0 or xy == 1);
    return d_number_exit[a][xy];
  }

  Point enter(u_int a, u_int i) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(i < d_exit[a].size());
    return d_enter[a][i];
  }

  Point exit(u_int a, u_int i) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(i < d_exit[a].size());
    return d_exit[a][i];
  }

  void display_tracks() const;

protected:

  /// \name Private Data
  /// \{

  /// Polar quadrature
  SP_polar d_polar;

  /// Number of azimuths per octant
  int d_number_azimuths_octant;

  /// Number of polar angles
  int d_number_polar;

  /// Track entrance points.
  vec2_point d_enter;

  /// Track exit points.
  vec2_point d_exit;

  /// Number of entrance points per angle per side.
  vec2_int d_number_enter;

  /// Number of exit points per angle per side.
  vec2_int d_number_exit;

  /// \}

};

} // end namespace detran

#endif /* QUADRATUREMOC_HH_ */

//---------------------------------------------------------------------------//
//              end of file QuadratureMOC.hh
//---------------------------------------------------------------------------//

