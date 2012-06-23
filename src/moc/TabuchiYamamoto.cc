/*
 * TabuchiYamamoto.cc
 *
 *  Created on: Jun 22, 2012
 *      Author: robertsj
 */

// Detran
#include "TabuchiYamamoto.hh"

namespace detran
{

TabuchiYamamoto::TabuchiYamamoto(int number_polar)
  : PolarQuadrature(number_polar)
{
  if (number_polar == 1)
  {
    d_sin[0]    = 0.798184;
    d_cos[0]    = 0.602413730042734;
    d_weight[0] = 1.0;
  }
  else if (number_polar == 2)
  {
    d_sin[0]    = 0.363900;
    d_cos[0]    = 0.931438023;
    d_sin[1]    = 0.899900;
    d_cos[1]    = 0.436096308;
    d_weight[0] = 0.212854;
    d_weight[1] = 0.787146;
  }
  else if (number_polar == 3)
  {
    d_sin[0]    = 0.166648;
    d_cos[0]    = 0.986016452;
    d_sin[1]    = 0.537707;
    d_cos[1]    = 0.84313177;
    d_sin[2]    = 0.932954;
    d_cos[2]    = 0.359995603;
    d_weight[0] = 0.046233;
    d_weight[1] = 0.283619;
    d_weight[2] = 0.670148;
  }
  else
  {
    THROW("Unsupported number of polar angles requested.");
  }

}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
