//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  QuadratureFactory.hh
 *  @brief QuadratureFactory class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_QUADRATUREFACTORY_HH_
#define detran_angle_QUADRATUREFACTORY_HH_

#include "BaseQuadrature.hh"
#include "Quadrature.hh"
#include "utilities/InputDB.hh"
#include <string>

namespace detran_angle
{

/**
 *  @class QuadratureFactory
 *  @brief Constructs quadratures.
 */
class ANGLE_EXPORT QuadratureFactory
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef BaseQuadrature::SP_basequadrature     SP_basequadrature;

  typedef detran_utilities::InputDB::SP_input   SP_input;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Build an angular quadrature set.
   *  @param input      User option database
   *  @param dimension  Quadrature dimension
   */
  static SP_quadrature build(SP_input input, const int dimension);

  /**
   *  @brief Build a 1-d base quadrature set
   *  @param    type  Name of quadrature type desired
   *  @param    n     Number of abscissa
   *  @param    a     Lower bound
   *  @param    b     Upper bound
   */
  static SP_basequadrature build_base(const std::string &type = "gl",
                                      const size_t       n   =  10,
                                      const double       a   = -1.0,
                                      const double       b   =  1.0);

};

} // end namespace detran_angle

#endif /* detran_angle_QUADRATUREFACTORY_HH_ */
