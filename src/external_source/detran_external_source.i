//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_external_source.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran external source.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") external_source
%{
#include <stddef.h>
#include "external_source/ExternalSource.hh"
#include "external_source/ConstantSource.hh"
#include "external_source/DiscreteSource.hh"
#include "external_source/IsotropicSource.hh"
%}

%feature("autodoc", "3");

// Hide templates from SWIG
%inline
{
#define EXTERNAL_SOURCE_EXPORT
#define EXTERNAL_SOURCE_TEMPLATE_EXPORT(...)
#define EXTERNAL_SOURCE_INSTANTIATE_EXPORT(...)
}

%include <pycontainer.swg>
%import "detran_utilities.i"
%import "detran_geometry.i"
%import "detran_angle.i"

%include "ExternalSource.hh"
%include "ConstantSource.hh"
%include "DiscreteSource.hh"
%include "IsotropicSource.hh"

%template(ExternalSourceSP)  detran_utilities::SP<detran_external_source::ExternalSource>;
%template(ConstantSourceSP)  detran_utilities::SP<detran_external_source::ConstantSource>;
%template(DiscreteSourceSP)  detran_utilities::SP<detran_external_source::DiscreteSource>;
%template(IsotropicSourceSP) detran_utilities::SP<detran_external_source::IsotropicSource>;

%template(vec_source) std::vector<detran_utilities::SP<detran_external_source::ExternalSource> >;

//---------------------------------------------------------------------------//
//              end of detran_external_source.i
//---------------------------------------------------------------------------//





