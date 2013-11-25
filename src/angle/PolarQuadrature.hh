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
// Various base quadratures
#include "angle/BaseGL.hh"
#include "angle/BaseGC.hh"
#include "angle/BaseUniform.hh"
#include "angle/TabuchiYamamoto.hh"
#include "angle/AbuShumaysDoubleRange.hh"
#include "angle/TriGauss.hh"

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

  /// Constructor, with optional normalization of weights to unity
  PolarQuadrature(const size_t number_polar, bool normalize = false);
  /// Number of polar angles per half space
  size_t number_polar() const;
  /// Get polar sines
  const vec_dbl& sin_theta() const;
  /// Get polar cosines
  const vec_dbl& cos_theta() const;
  /// Return base name
  static std::string name() {return B::name();}

private:

  /// Number of polar angles for a single half space.
  size_t d_number_polar;
  /// Vector of polar cosines
  vec_dbl d_cos_theta;
  /// Vector of polar sines
  vec_dbl d_sin_theta;

};

//----------------------------------------------------------------------------//
// CONVENIENCE TYPEDEFS
//----------------------------------------------------------------------------//

typedef PolarQuadrature<BaseGL>                 PolarGL;
typedef PolarQuadrature<BaseDGL>                PolarDGL;
typedef PolarQuadrature<BaseGC>                 PolarGC;
typedef PolarQuadrature<BaseDGC>                PolarDGC;
typedef PolarQuadrature<BaseUniform>            PolarU;
typedef PolarQuadrature<BaseUniformCosine>      PolarUC;
typedef PolarQuadrature<BaseSimpson>            PolarS;
typedef PolarQuadrature<TabuchiYamamoto>        PolarTY;
typedef PolarQuadrature<AbuShumaysDoubleRange>  PolarASDR;
typedef PolarQuadrature<TriGauss>               PolarTG;

} // end namespace detran_angle

#endif /* detran_angle_POLARQUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of PolarQuadrature.hh
//----------------------------------------------------------------------------//
