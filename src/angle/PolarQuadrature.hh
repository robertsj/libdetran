//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PolarQuadrature.hh
 *  @brief PolarQuadrature class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_POLARQUADRATURE_HH_
#define detran_angle_POLARQUADRATURE_HH_

#include "angle/Quadrature.hh"

namespace detran_angle
{

/**
 *  @class PolarQuadrature
 *  @brief Defines a quadrature over [-1, 1]
 *
 *  Using any 1-D quadrature implementing @ref BaseQuadrature, this class
 *  yields a complete angular quadrature for slab geometry or for use in
 *  defining a @ref ProductQuadrature.
 *
 *  @tparam B   Base quadrature class
 */
template <class B>
class ANGLE_EXPORT PolarQuadrature: public Quadrature
{

public:

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    number_polar   Number of angles in positive half space
   *  @param    normalize      Flag to normalize quadrant weights to unity
   */
  PolarQuadrature(const size_t number_polar, bool normalize = false);

  /// Number of polar angles per half space
  size_t number_polar() const;

  /**
   *  @brief Return a polar sine
   *  @param    p  Polar index in octant
   */
  double sin_theta(const size_t p) const;

  /**
   *  @brief Return a polar cosine
   *  @param    p  Polar index in octant
   */
  double cos_theta(const size_t p) const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Number of polar angles for a single half space.
  size_t d_number_polar;
  /// Vector of polar sines
  std::vector<double> d_sin_theta;
  /// Vector of polar cosines
  std::vector<double> d_cos_theta;

};

} // end namespace detran_angle

#endif /* POLARQUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of PolarQuadrature.hh
//----------------------------------------------------------------------------//
