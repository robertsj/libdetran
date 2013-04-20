//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_angle.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran angular system.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") angle
%{
#include <stddef.h>
#include "angle/detran_angle.hh"
%}

%import "utilities/detran_utilities.i"

// Base angle classes and utilities
%include "Quadrature.hh"
%include "ProductQuadrature.hh"
%include "QuadratureMOC.hh"
%include "MomentToDiscrete.hh"
%include "QuadratureFactory.hh"
%include "MomentIndexer.hh"

%template(QuadratureSP)         detran_utilities::SP<detran_angle::Quadrature>;
%template(ProductQuadratureSP)  detran_utilities::SP<detran_angle::ProductQuadrature>;
%template(QuadratureMOCSP)      detran_utilities::SP<detran_angle::QuadratureMOC>;
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

// Cast from base quadrature to moc quadrature
detran_utilities::SP<detran_angle::QuadratureMOC> 
cast_base_to_moc(detran_utilities::SP<detran_angle::Quadrature>* q)
{
  return detran_utilities::SP<detran_angle::QuadratureMOC>(*q);
} 

} // end inline



//---------------------------------------------------------------------------//
//              end of detran_angle.i
//---------------------------------------------------------------------------//
