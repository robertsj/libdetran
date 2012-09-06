//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   QuadratureFactory.hh
 * \author robertsj
 * \date   Apr 11, 2012
 * \brief  QuadratureFactory class definition.
 */
//---------------------------------------------------------------------------//

#ifndef QUADRATUREFACTORY_HH_
#define QUADRATUREFACTORY_HH_

#include "Quadrature.hh"
#include "utilities/InputDB.hh"
#include <string>

namespace detran_angle
{

/*!
 *  \class QuadratureFactory
 *  \brief Constructs quadratures.
 */
class QuadratureFactory
{

public:

  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef detran_utilities::InputDB::SP_input   SP_input;

  /*!
   *  \brief Build a quadrature set.
   *
   *  \param q          Smart pointer to quadrature being constructed
   *  \param input      Smart point to input.
   *  \param dimension  Quadrature dimension
   */
  void build(SP_quadrature &q, SP_input input, int dimension);


  /*!
   *  \brief Build a quadrature set.
   *
   *  \param q          Smart pointer to quadrature being constructed
   *  \param type       Quadrature type
   *  \param order      Quadrature order
   *  \param dimension  Quadrature dimension
   */
  void build(SP_quadrature &q, std::string type, int order, int dimension);
//
//  /*!
//   *  \brief Build an MOC quadrature set.
//   *
//   *  \param q                          Smart pointer to quadrature being constructed
//   *  \param type                       Azimuthal quadrature type
//   *  \param polar                      Polar quadrature type
//   *  \param number_azimuths_octant     Number of azimuths per octant
//   *  \param number_polar_octant        Number of polar angles per octant
//   */
//  void build_MOC(SP_quadrature &q,
//             std::string type,
//             std::string polar,
//             int number_azimuths_octant,
//             int number_polar_octant);

};

} // end namespace detran_angle

#endif /* QUADRATUREFACTORY_HH_ */
