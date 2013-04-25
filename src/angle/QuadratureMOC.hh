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

#include "angle/angle_export.hh"
#include "Quadrature.hh"
#include "PolarQuadrature.hh"
#include "utilities/Point.hh"
#include <cmath>

namespace detran_angle
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
 *  \todo It may be useful to define an interface that would
 *        require the polar and azimith components be
 *        constructed separately with an inherited method
 *        to build the product set.  This would enforce a
 *        well-quantified ordering for future analysis.
 *
 */
class ANGLE_EXPORT QuadratureMOC: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<QuadratureMOC>           SP_quadrature;
  typedef Quadrature                                    Base;
  typedef Base::SP_quadrature                           SP_base;
  typedef detran_utilities::SP<PolarQuadrature>         SP_polar;
  typedef detran_utilities::Point                       Point;
  typedef std::vector<Point>                            vec_point;
  typedef std::vector<vec_point>                        vec2_point;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec2_int                    vec2_int;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param dim                  Problem dimension (2 or 3d)
   *  \param num_azimuths_octant  Number azimuths per octant
   *  \param num_polar            Number of polar angles
   *  \param name                 Name of this azimuthal quadrature
   *  \param polar                Name of the requested polar quadrature
   */
  QuadratureMOC(size_t dim,
                size_t num_azimuths_octant,
                size_t num_polar,
                std::string name,
                std::string polar);

  size_t number_azimuths_octant() const
  {
    return d_number_azimuths_octant;
  }

  size_t number_polar_octant() const
  {
    return d_number_polar;
  }

  size_t number_tracks(size_t a)
  {
    Require(a < 2*d_number_azimuths_octant);
    return d_exit[a].size();
  }

  size_t number_enter(size_t a, size_t xy) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(xy == 0 or xy == 1);
    return d_number_enter[a][xy];
  }

  int number_exit(size_t a, size_t xy) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(xy == 0 or xy == 1);
    return d_number_exit[a][xy];
  }

  Point enter(size_t a, size_t i) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(i < d_exit[a].size());
    return d_enter[a][i];
  }

  Point exit(size_t a, size_t i) const
  {
    Require(a < 2*d_number_azimuths_octant);
    Require(i < d_exit[a].size());
    return d_exit[a][i];
  }

  /// Return angle within octant given azimuth and polar.
  size_t angle(size_t a, size_t p) const
  {
    Require(a < d_number_azimuths_octant);
    Require(p < d_number_polar);
    return p + a * d_number_polar;
  }

  /// Return azimuth index for angle within octant
  size_t azimuth(size_t angle) const
  {
    Require(angle < d_number_angles_octant);
    double tmp = double(angle % d_number_angles_octant)/double(d_number_polar);
    return int(std::floor(tmp));
  }

  /// Return polar index for cardinal index
  size_t polar(size_t angle) const
  {
    Require(angle < d_number_angles_octant);
    return angle % d_number_polar;
  }

  double sin_theta(size_t p) const
  {
    Require(p < d_number_polar);
    return d_polar->sin_theta(p);
  }

  double cos_theta(size_t p) const
  {
    Require(p < d_number_polar);
    return d_polar->cos_theta(p);
  }

  double phi(size_t a) const
  {
    Require(a < 2 * d_number_azimuths_octant);
    return d_phi[a];
  }

  double sin_phi(size_t a) const
  {
    Require(a < 2 * d_number_azimuths_octant);
    return d_sin_phi[a];
  }

  double cos_phi(size_t a) const
  {
    Require(a < 2 * d_number_azimuths_octant);
    return d_cos_phi[a];
  }

  double spacing(size_t a) const
  {
    Require(a < 2 * d_number_azimuths_octant);
    return d_spacing[a];
  }

  double azimuth_weight(size_t a) const
  {
    Require(a < 2 * d_number_azimuths_octant);
    return d_azimuth_weight[a];
  }

  void display_tracks() const;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Polar quadrature
  SP_polar d_polar;
  /// Azimuth angles
  vec_dbl d_phi;
  /// Azimuth cosines
  vec_dbl d_cos_phi;
  /// Azimuth sines
  vec_dbl d_sin_phi;
  /// Normalized spacing (one 1x1 box)
  vec_dbl d_spacing;
  /// Azimuth weights
  vec_dbl d_azimuth_weight;
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

};

ANGLE_TEMPLATE_EXPORT(detran_utilities::SP<QuadratureMOC>)

} // end namespace detran

#endif /* QUADRATUREMOC_HH_ */

//---------------------------------------------------------------------------//
//              end of file QuadratureMOC.hh
//---------------------------------------------------------------------------//

