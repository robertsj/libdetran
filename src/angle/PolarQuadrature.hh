//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolarQuadrature.hh
 * \brief  PolarQuadrature class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef POLARQUADRATURE_HH_
#define POLARQUADRATURE_HH_

#include "angle/angle_export.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_angle
{

/*!
 *  \class PolarQuadrature
 *  \brief Polar quadrature definition, nominally for MOC
 */
class ANGLE_EXPORT PolarQuadrature
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<PolarQuadrature>       SP_polar;
  typedef detran_utilities::size_t                    size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param number_polar   Number of angles in positive half space
   */
  PolarQuadrature(size_t number_polar)
    : d_number_polar(number_polar)
    , d_sin(number_polar, 0.0)
    , d_cos(number_polar, 0.0)
    , d_weight(number_polar, 0.0)
  {
    Require(number_polar > 0);
  }

  /// Virtual destructor.
  virtual ~PolarQuadrature(){};

  size_t number_polar() const
  {
    return d_number_polar;
  }

  /*!
   *  \brief Return a polar sine
   *  \param p  Polar index in octant
   */
  double sin_theta(size_t p) const
  {
    Require(p < d_number_polar);
    return d_sin[p];
  }

  /*!
   *  \brief Return a polar cosine
   *  \param p  Polar index in octant
   */
  double cos_theta(size_t p) const
  {
    Require(p < d_number_polar);
    return d_cos[p];
  }

  /*!
   *  \brief Return a polar weight
   *  \param p  Polar index in octant
   */
  double weight(size_t p) const
  {
    Require(p < d_number_polar);
    return d_weight[p];
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Number of polar angles for a single half space.
  size_t d_number_polar;
  /// Vector of sines
  std::vector<double> d_sin;
  /// Vector of cosines
  std::vector<double> d_cos;
  /// Vector of weights (sums to 1).
  std::vector<double> d_weight;

};

ANGLE_TEMPLATE_EXPORT(detran_utilities::SP<PolarQuadrature>)

} // end namespace detran_angle

#endif /* POLARQUADRATURE_HH_ */
