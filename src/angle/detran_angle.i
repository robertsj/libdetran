//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_angle.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran angular system.
 */
//---------------------------------------------------------------------------//

%include "utilities/detran_utilities.i"

// Base angle classes
%include "Quadrature.hh"
//
%include "GaussLegendre.hh"
%include "GaussChebyshev.hh"
%include "DPN.hh"
%include "DTN.hh"
//
%include "LevelSymmetric.hh"
%include "QuadrupleRange.hh"
%include "UniformEqual.hh"
%include "MomentToDiscrete.hh"
// MOC angle classes
%include "PolarQuadrature.hh"
%include "TabuchiYamamoto.hh"
%include "QuadratureMOC.hh"
%include "Collocated.hh"
%include "Uniform.hh"
// Product quadrature
%include "ProductQuadrature.hh"
%include "ChebyshevLegendre.hh"

%template(QuadratureSP)        detran_utilities::SP<detran_angle::Quadrature>;
%template(ProductQuadratureSP) detran_utilities::SP<detran_angle::ProductQuadrature>;
//
%template(GaussLegendreSP)    detran_utilities::SP<detran_angle::GaussLegendre>;
%template(GaussChebyshevSP)   detran_utilities::SP<detran_angle::GaussChebyshev>;
%template(DPN_SP)             detran_utilities::SP<detran_angle::DPN>;
%template(DTN_SP)             detran_utilities::SP<detran_angle::DTN>;
//
%template(LevelSymmetricSP)   detran_utilities::SP<detran_angle::LevelSymmetric>;
%template(QuadrupleRangeSP)   detran_utilities::SP<detran_angle::QuadrupleRange>;
%template(UniformEqualSP)     detran_utilities::SP<detran_angle::UniformEqual>;
//
%template(MomentToDiscreteSP) detran_utilities::SP<detran_angle::MomentToDiscrete>;

%template(QuadratureMOCSP)    detran_utilities::SP<detran_angle::QuadratureMOC>;
%template(CollocatedSP)       detran_utilities::SP<detran_angle::Collocated>;
%template(UniformSP)          detran_utilities::SP<detran_angle::Uniform>;
//
%template(PolarQuadratureSP)  detran_utilities::SP<detran_angle::PolarQuadrature>;
%template(TabuchiYamamotoSP)  detran_utilities::SP<detran_angle::TabuchiYamamoto>;

//%inline
//{
//
//// Cast from base quadrature to product quadrature
//detran_utilities::SP<detran_angle::ProductQuadrature> 
//cast_base_to_product(detran_utilities::SP<detran_angle::Quadrature> >* q)
//{
//  return detran_utilities::SP<detran_angle::ProductQuadrature>(*q);
//} 
//
//}

//---------------------------------------------------------------------------//
//              end of detran_angle.i
//---------------------------------------------------------------------------//
