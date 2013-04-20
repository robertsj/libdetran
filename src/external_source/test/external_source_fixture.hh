//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   external_source_fixture.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  external_source_fixture class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Setup
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"

// Detran
#include "ConstantSource.hh"
#include "DiscreteSource.hh"
#include "IsotropicSource.hh"

namespace detran_test
{
  typedef detran_external_source::
    ExternalSource::SP_externalsource SP_externalsource;

  static SP_externalsource constant_source_fixture(int dimension)
  {

    SP_mesh mesh;
    SP_quadrature quad;
    if (dimension == 1)
    {
      //
    }
    else if(dimension == 2)
    {
      mesh = mesh_2d_fixture();
      quad = quadruplerange_fixture();
    }
    else
    {
      mesh = mesh_3d_fixture();
    }

    SP_externalsource
      q_e(new detran_external_source::ConstantSource(2, mesh, 1.0, quad));
    return q_e;

  }


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of external_source_fixture.hh
//---------------------------------------------------------------------------//
