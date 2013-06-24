//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TabuchiYamamoto.hh
 *  @brief TabuchiYamamoto class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_TABUCHIYAMAMOTO_HH_
#define detran_angle_TABUCHIYAMAMOTO_HH_

#include "angle/BaseQuadrature.hh"

namespace detran_angle
{

/**
 *  @class TabuchiYamamoto
 *  @brief Quadrature of Tabuchi and Yamamoto for polar integration
 *
 *  Tabuchi and Yamamoto define an optimal polar quadrature such that
 *  the maximum error of an approximate Ki(3,x) over all x is minimized.
 *  Ki(3,x) is a Bickley-Naylor function, which is related to modified
 *  Bessel functions.
 *
 *  Their original paper gives quadratures for 1, 2, and 3 abscissa.
 *  We extend that to 5 abscissa.  Additionally, we implement another
 *  size abscissa based on an L2 minimization of the error (rather
 *  than the infinity norm)... not yet implemented.
 *
 */
class ANGLE_EXPORT TabuchiYamamoto: public BaseQuadrature
{

public:

  /// Return the quadrature name.
  static std::string name() {return "ty";}

private:

  void build_impl(c_dbl a, c_dbl b);

};

} // end namespace detran

#endif /* detran_angle_TABUCHIYAMAMOTO_HH_ */

//----------------------------------------------------------------------------//
//              end of file TabuchiYamamoto.hh
//----------------------------------------------------------------------------//
