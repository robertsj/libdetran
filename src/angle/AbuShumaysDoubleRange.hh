//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  AbuShumaysDoubleRange.hh
 *  @brief AbuShumaysDoubleRange class definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_ABUSHUMAYSDOUBLERANGE_HH_
#define detran_angle_ABUSHUMAYSDOUBLERANGE_HH_

#include "angle/BaseQuadrature.hh"

namespace detran_angle
{

class ANGLE_EXPORT AbuShumaysDoubleRange: public BaseQuadrature
{

public:

  static std::string name() {return "asdr";}

private:

  void build_impl(c_dbl a, c_dbl b);

  /// Get the *cosine* of the polar angle or its weight
  double get_theta(const size_t N, const size_t i, const size_t j);

};

} // end namespace detran_angle

#endif /* detran_angle_ABUSHUMAYSDOUBLERANGE_HH_ */

//----------------------------------------------------------------------------//
//              end of AbuShumaysDoubleRange.hh
//----------------------------------------------------------------------------//
