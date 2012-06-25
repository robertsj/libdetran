//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadraturemoc_fixture.hh
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Quadrature fixtures for testing.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//


#ifndef QUADRATUREMOC_FIXTURE_HH_
#define QUADRATUREMOC_FIXTURE_HH_

// Detran includes
#include "QuadratureMOC.hh"
#include "Unfiorm.hh"

// Detran utilities
#include "DBC.hh"

// System includes

namespace detran_test
{

typedef detran::QuadratureMOC::SP_quadrature SP_quadrature;


/*!
 *  \brief Create a Uniform MOC quadrature for testing.
 *
 *  Note, for classes like this, it's nearly as easy simply to
 *  create directly when required.  However, setting it once
 *  here ensures consistency at all use points.
 */
static SP_quadrature uniform_fixture()
{
  // Define quadrature.
  SP_quadrature q;
  q = new detran::QuadrupleRange(2, 2);

  // Return the fixture.
  return q;
}

/*!
 *  \brief Create a GaussLegendre quadrature for testing.
 *
 *  Note, for classes like this, it's nearly as easy simply to
 *  create directly when required.  However, setting it once
 *  here ensures consistency at all use points.
 */
static SP_quadrature gausslegendre_fixture()
{
  // Define quadrature.
  SP_quadrature q;
  q = new detran::GaussLegendre(8);

  // Return the fixture.
  return q;
}


} // end namespace detran_test

#endif /* QUADRATUREMOC_FIXTURE_HH_ */

//---------------------------------------------------------------------------//
//              end of quadraturemoc_fixture.hh
//---------------------------------------------------------------------------//
