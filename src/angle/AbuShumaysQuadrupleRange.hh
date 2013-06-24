//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  AbuShumaysQuadrupleRange.hh
 *  @brief AbuShumaysQuadrupleRange
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_ABUSHUMAYSQUADRUPLERANGE_HH_
#define detran_angle_ABUSHUMAYSQUADRUPLERANGE_HH_

#include "angle/BaseQuadrature.hh"

namespace detran_angle
{

class ANGLE_EXPORT AbuShumaysQuadrupleRange: public BaseQuadrature
{

public:

  static std::string name() {return "asqr";}

private:

  void build_impl(c_dbl a, c_dbl b);

  /// Get the azimuth or its weight
  double get_phi(const size_t N, const size_t i, const size_t j);

};

} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of AbuShumaysQuadrupleRange.hh
//----------------------------------------------------------------------------//

#endif /* detran_angle_ABUSHUMAYSQUADRUPLERANGE_HH_ */
