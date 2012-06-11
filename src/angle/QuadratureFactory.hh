//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   QuadratureFactory.hh
 * \author robertsj
 * \date   Apr 11, 2012
 * \brief  QuadratureFactory class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef QUADRATUREFACTORY_HH_
#define QUADRATUREFACTORY_HH_

// Detran
#include "Quadrature.hh"

// System
#include <string>

namespace detran
{

/*!
 *  \class QuadratureFactory
 *  \brief Constructs quadratures.
 */
class QuadratureFactory
{

public:

  /// Quadrature type.
  typedef Quadrature::SP_quadrature SP_quadrature;

  /*!
   *  \brief Build a quadrature set.
   *
   *  \param q          Smart pointer to quadrature being constructed
   *  \param type       Quadrature type
   *  \param order      Quadrature order
   *  \param dimension  Quadrature dimension
   */
  void build(SP_quadrature &q, std::string type, int order, int dimension);
};

}

#endif /* QUADRATUREFACTORY_HH_ */
