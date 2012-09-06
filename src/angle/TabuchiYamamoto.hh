//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TabuchiYamamoto.hh
 * \brief  TabuchiYamamoto class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef TABUCHIYAMAMOTO_HH_
#define TABUCHIYAMAMOTO_HH_

#include "PolarQuadrature.hh"

namespace detran_angle
{

/*!
 *  \class TabuchiYamamoto
 *  \brief Quadrature of Tabuchi and Yamamoto for polar integration
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
