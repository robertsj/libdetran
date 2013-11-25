//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  ConstantSource.cc
 *  @brief ConstantSource member definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//---------------------------------------------------------------------------//

#include "ConstantSource.hh"

namespace detran_external_source
{

//---------------------------------------------------------------------------//
ConstantSource::ConstantSource(size_t number_groups,
                               SP_mesh mesh,
                               double source,
                               SP_quadrature quadrature)
  : ExternalSource(number_groups, mesh, quadrature)
  , d_source(source)
  , d_discrete_source(source *
      detran_angle::Quadrature::angular_norm(mesh->dimension()))
{
  /* ... */
}

//---------------------------------------------------------------------------//
ConstantSource::SP_externalsource
ConstantSource::Create(size_t         number_groups,
                       SP_mesh        mesh,
                       double         source,
                       SP_quadrature  quadrature)
{
  SP_externalsource
    p(new ConstantSource(number_groups, mesh, source, quadrature));
  return p;
}

} // end namespace detran_external_source

//---------------------------------------------------------------------------//
//              end of file ConstantSource.cc
//---------------------------------------------------------------------------//
