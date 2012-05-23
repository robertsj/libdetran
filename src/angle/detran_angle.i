//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_angle.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran angular system.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "Quadrature.hh"
%include "GaussLegendre.hh"
%include "LevelSymmetric.hh"
%include "QuadrupleRange.hh"
%include "UniformEqual.hh"
%include "MomentToDiscrete.hh"

%template(QuadratureSP)     detran::SP<detran::Quadrature>;
%template(GaussLegendreSP)  detran::SP<detran::GaussLegendre>;
%template(LevelSymmetricSP) detran::SP<detran::LevelSymmetric>;
%template(QuadrupleRangeSP) detran::SP<detran::QuadrupleRange>;
%template(UniformEqualSP)   detran::SP<detran::UniformEqual>;
//
%template(MomentToDiscrete1D) detran::MomentToDiscrete<detran::_1D>;
%template(MomentToDiscrete2D) detran::MomentToDiscrete<detran::_2D>;
%template(MomentToDiscrete3D) detran::MomentToDiscrete<detran::_3D>;
//
%template(MomentToDiscrete1DSP) detran::SP<detran::MomentToDiscrete<detran::_1D> >;
%template(MomentToDiscrete2DSP) detran::SP<detran::MomentToDiscrete<detran::_2D> >;
%template(MomentToDiscrete3DSP) detran::SP<detran::MomentToDiscrete<detran::_3D> >;


//---------------------------------------------------------------------------//
//              end of detran_angle.i
//---------------------------------------------------------------------------//
