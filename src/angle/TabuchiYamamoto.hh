//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TabuchiYamamoto.hh
 *  @brief  TabuchiYamamoto class definition.
 *  @author Jeremy Roberts
 *  @date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef TABUCHIYAMAMOTO_HH_
#define TABUCHIYAMAMOTO_HH_

#include "PolarQuadrature.hh"

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
 *  We extend that to 6 abscissa.  Additionally, we implement another
 *  size abscissa based on an L2 minimization of the error (rather
 *  than the infinity norm).
 *
 */
class TabuchiYamamoto: public PolarQuadrature
{

public:

  explicit TabuchiYamamoto(size_t number_polar);

  /// SP Constructor.
  static detran_utilities::SP<PolarQuadrature>
  Create(size_t num_polar)
  {
    detran_utilities::SP<TabuchiYamamoto>
      p(new TabuchiYamamoto(num_polar));
    return p;
  }

  ~TabuchiYamamoto(){};

private:

};

} // end namespace detran

#endif /* TABUCHIYAMAMOTO_HH_ */

//---------------------------------------------------------------------------//
//              end of file TabuchiYamamoto.hh
//---------------------------------------------------------------------------//
