//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature_fixture.hh
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Quadrature fixtures for testing.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef QUADRATURE_FIXTURE_HH_
#define QUADRATURE_FIXTURE_HH_

// Detran includes
#include "QuadrupleRange.hh"
#include "GaussLegendre.hh"

// Detran utilities
#include "DBC.hh"

// System includes

namespace detran_test
{

typedef detran::Quadrature::SP_quadrature SP_quadrature;


/*!
 *  \brief Create a QuadrupleRange quadrature for testing.
 *
 *  Note, for classes like this, it's nearly as easy simply to
 *  create directly when required.  However, setting it once
 *  here ensures consistency at all use points.
 */
static SP_quadrature quadruplerange_fixture()
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

#endif /* QUADRATURE_FIXTURE_HH_ */

//---------------------------------------------------------------------------//
//              end of quadrature_fixture.hh
//---------------------------------------------------------------------------//
