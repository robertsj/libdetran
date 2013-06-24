//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  AzimuthalQuadrature.hh
 *  @brief AzimuthalQuadrature class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_AZIMUTHALQUADRATURE_HH_
#define detran_angle_AZIMUTHALQUADRATURE_HH_

#include "angle/Quadrature.hh"
// Various base quadratures
#include "angle/BaseGL.hh"
#include "angle/BaseUniform.hh"
#include "angle/AbuShumaysQuadrupleRange.hh"

namespace detran_angle
{

/**
 *  @class AzimuthalQuadrature
 *  @brief Defines a quadrature over [0, pi/2]
 *
 *  Unlike @ref PolarQuadrature, this class is not a complete transport
 *  quadrature.  It is only used for defining product quadratures.
 *
 *  @tparam B   Base quadrature class
 */
template <class B>
class ANGLE_EXPORT AzimuthalQuadrature
{

public:

  typedef Quadrature::vec_dbl   vec_dbl;
  typedef Quadrature::size_t    size_t;

  /// Constructor, with optional normalization of weights to @f$ \pi/2 @f$
  AzimuthalQuadrature(const size_t number_azimuth, bool normalize = false);
  /// Number of polar angles per half space
  size_t number_azimuth() const;
  /// Get azimuths
  const vec_dbl& phi() const;
  /// Get azimuth cosines
  const vec_dbl& cos_phi() const;
  /// Get azimuth sines
  const vec_dbl& sin_phi() const;
  /// Get all azimuthal weights
  const vec_dbl& weights() const;
  /// Return base name
  static std::string name() {return B::name();}

private:

  /// Number of azimuths per quadrant
  size_t d_number_azimuth;
  /// Vector of azimuths
  vec_dbl d_phi;
  /// Vector of azimuths
  vec_dbl d_cos_phi;
  /// Vector of azimuths
  vec_dbl d_sin_phi;
  /// Vector of weights
  vec_dbl d_weight;

};

//----------------------------------------------------------------------------//
// CONVENIENCE TYPEDEFS
//----------------------------------------------------------------------------//

typedef AzimuthalQuadrature<BaseGL>                     AzimuthalGL;
typedef AzimuthalQuadrature<BaseDGL>                    AzimuthalDGL;
typedef AzimuthalQuadrature<BaseUniform>                AzimuthalU;
typedef AzimuthalQuadrature<BaseSimpson>                AzimuthalS;
typedef AzimuthalQuadrature<AbuShumaysQuadrupleRange>   AzimuthalASQR;

} // end namespace detran_angle

#endif /* detran_angle_AZIMUTHALQUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of AzimuthalQuadrature.hh
//----------------------------------------------------------------------------//
