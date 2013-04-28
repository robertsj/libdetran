//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  TabuchiYamamoto.cc
 *  @brief TabuchiYamamoto member definitions
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//---------------------------------------------------------------------------//

#include "TabuchiYamamoto.hh"
#include <cmath>

namespace detran_angle
{

//---------------------------------------------------------------------------//
TabuchiYamamoto::TabuchiYamamoto(size_t number_polar)
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
  else if (number_polar == 4)
  {
    d_sin[0] = 0.08370291;
    d_sin[1] = 0.29149397;
    d_sin[2] = 0.63350448;
    d_sin[3] = 0.94935787;
    for (int i = 0; i < 4; ++i)
      d_cos[i] = std::sqrt(1.0 - d_sin[i]*d_sin[i]);
    d_weight[0] = 0.01173545;
    d_weight[1] = 0.08582263;
    d_weight[2] = 0.30817276;
    d_weight[3] = 0.59426914;

  }
  else if (number_polar == 5)
  {
    d_sin[0] = 0.04682256;
    d_sin[1] = 0.16656767;
    d_sin[2] = 0.39033439;
    d_sin[3] = 0.70255004;
    d_sin[4] = 0.96043925;
    for (int i = 0; i < 5; ++i)
      d_cos[i] = std::sqrt(1.0 - d_sin[i]*d_sin[i]);
    d_weight[0] = 0.00368232000000002;
    d_weight[1] = 0.02830905000000000;
    d_weight[2] = 0.11831753000000000;
    d_weight[3] = 0.31699287999999998;
    d_weight[4] = 0.53269822000000000;
  }
  else
  {
    THROW("Unsupported number of polar angles requested.");
  }

}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of file TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
