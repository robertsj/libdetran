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
  typedef detran::ExternalSource::SP_source SP_source;

  static SP_source constant_source_fixture(int dimension)
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

    SP_source q_e;
    q_e = new detran::ConstantSource(mesh, quad, 2);
    return q_e;

  }


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of external_source_fixture.hh
//---------------------------------------------------------------------------//
