//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   detran_angle.i
 *  @brief  Python interface for detran angular system
 *  @note   Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") angle
%{
#include <stddef.h>
#include "angle/detran_angle.hh"
%}

%feature("autodoc", "3");

// Hide templates from SWIG
%inline
{
#define ANGLE_EXPORT
#define ANGLE_TEMPLATE_EXPORT(...)
#define ANGLE_INSTANTIATE_EXPORT(...)
}

%import "utilities/detran_utilities.i"

// Base angle classes and utilities
%include "BaseQuadrature.hh"
%include "Quadrature.hh"
%include "ProductQuadrature.hh"
%include "MomentToDiscrete.hh"
%include "QuadratureFactory.hh"
%include "MomentIndexer.hh"

%template(BaseQuadratureSP)     detran_utilities::SP<detran_angle::BaseQuadrature>;
%template(QuadratureSP)         detran_utilities::SP<detran_angle::Quadrature>;
%template(ProductQuadratureSP)  detran_utilities::SP<detran_angle::ProductQuadrature>;
%template(MomentToDiscreteSP)   detran_utilities::SP<detran_angle::MomentToDiscrete>;
%template(MomentIndexerSP)      detran_utilities::SP<detran_angle::MomentIndexer>;

%inline
{

// Cast from base quadrature to product quadrature
detran_utilities::SP<detran_angle::ProductQuadrature> 
cast_base_to_product(detran_utilities::SP<detran_angle::Quadrature>* q)
{
  return detran_utilities::SP<detran_angle::ProductQuadrature>(*q);
} 

} // end inline

//----------------------------------------------------------------------------//
//              end of detran_angle.i
//----------------------------------------------------------------------------//
