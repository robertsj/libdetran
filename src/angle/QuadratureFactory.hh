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

#include "Quadrature.hh"

#include <string>

namespace detran
{

class QuadratureFactory
{

  public:

    //! Quadrature type.
    typedef Quadrature::SP_quadrature SP_quadrature;

    // Build a quadrature set.
    void build(SP_quadrature &q, std::string type, int order, int dimension);
};

}

#endif /* QUADRATUREFACTORY_HH_ */
