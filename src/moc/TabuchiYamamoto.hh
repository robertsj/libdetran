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

class TabuchiYamamoto: public PolarQuadrature
{

public:

  TabuchiYamamoto(int number_polar);

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
