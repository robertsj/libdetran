/*
 * TabuchiYamamoto.hh
 *
 *  Created on: Jun 22, 2012
 *      Author: robertsj
 */

#ifndef TABUCHIYAMAMOTO_HH_
#define TABUCHIYAMAMOTO_HH_

// Detran
#include "PolarQuadrature.hh"

namespace detran
{

/*!
 *  \class TabuchiYamamoto
 *  \brief Quadrature of Tabuchi and Yamamoto for polar integration
 */
class TabuchiYamamoto: public PolarQuadrature
{

public:

  explicit TabuchiYamamoto(int number_polar);

  /// SP Constructor.
  static SP<PolarQuadrature> Create(int num_polar)
  {
    SP<TabuchiYamamoto> p(new TabuchiYamamoto(num_polar));
    return p;
  }

  ~TabuchiYamamoto(){};

  bool is_valid() const
  {
    Require(d_number_polar > 0.0);
    Require(d_sin.size() == d_number_polar);
    return true;
  }

private:

};

} // end namespace detran

#endif /* TABUCHIYAMAMOTO_HH_ */

//---------------------------------------------------------------------------//
//              end of file TabuchiYamamoto.hh
//---------------------------------------------------------------------------//
