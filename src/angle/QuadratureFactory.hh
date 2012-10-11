//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   QuadratureFactory.hh
 *  @author robertsj
 *  @date   Apr 11, 2012
 *  @brief  QuadratureFactory class definition.
 */
//---------------------------------------------------------------------------//

#ifndef QUADRATUREFACTORY_HH_
#define QUADRATUREFACTORY_HH_

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
   *  \brief Build a quadrature set.
   *
   *  \param q          Smart pointer to quadrature being constructed
   *  \param input      Smart point to input.
   *  \param dimension  Quadrature dimension
   */
  void build(SP_quadrature &q, SP_input input, int dimension);

};

} // end namespace detran_angle

#endif /* QUADRATUREFACTORY_HH_ */
