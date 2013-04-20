//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   QuadratureFactory.hh
 *  @author robertsj
 *  @date   Apr 11, 2012
 *  @brief  QuadratureFactory class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_QUADRATUREFACTORY_HH_
#define detran_angle_QUADRATUREFACTORY_HH_

#include "Quadrature.hh"
#include "utilities/InputDB.hh"
#include <string>

namespace detran_angle
{

/**
 *  @class QuadratureFactory
 *  @brief Constructs quadratures.
 */
class QuadratureFactory
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef detran_utilities::InputDB::SP_input   SP_input;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Build a quadrature set.
   *
   *  @param q          Smart pointer to quadrature being constructed
   *  @param input      Smart point to input.
   *  @param dimension  Quadrature dimension
   */
  void build(SP_quadrature &q, SP_input input, const int dimension);

  //
  SP_quadrature build(SP_input input, const int dimension);

  /// Print out the available quadratures, etc.
  void help() const;

};

} // end namespace detran_angle

#endif /* detran_angle_QUADRATUREFACTORY_HH_ */
