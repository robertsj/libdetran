//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_utilities.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran angular system.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

// Base angle classes
%include "Quadrature.hh"
%include "GaussLegendre.hh"
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

%template(QuadratureSP)     detran_utilities::SP<detran_angle::Quadrature>;
%template(GaussLegendreSP)  detran_utilities::SP<detran_angle::GaussLegendre>;
%template(LevelSymmetricSP) detran_utilities::SP<detran_angle::LevelSymmetric>;
%template(QuadrupleRangeSP) detran_utilities::SP<detran_angle::QuadrupleRange>;
%template(UniformEqualSP)   detran_utilities::SP<detran_angle::UniformEqual>;
//
%template(MomentToDiscrete1D) detran_angle::MomentToDiscrete<detran::_1D>;
%template(MomentToDiscrete2D) detran_angle::MomentToDiscrete<detran::_2D>;
%template(MomentToDiscrete3D) detran_angle::MomentToDiscrete<detran::_3D>;
//
%template(MomentToDiscrete1DSP) detran_utilities::SP<detran_angle::MomentToDiscrete<detran::_1D> >;
%template(MomentToDiscrete2DSP) detran_utilities::SP<detran_angle::MomentToDiscrete<detran::_2D> >;
%template(MomentToDiscrete3DSP) detran_utilities::SP<detran_angle::MomentToDiscrete<detran::_3D> >;

%template(QuadratureMOCSP)    detran_utilities::SP<detran_angle::QuadratureMOC>;
%template(CollocatedSP)       detran_utilities::SP<detran_angle::Collocated>;
%template(UniformSP)          detran_utilities::SP<detran_angle::Uniform>;
//
%template(PolarQuadratureSP)  detran_utilities::SP<detran_angle::PolarQuadrature>;
%template(TabuchiYamamotoSP)  detran_utilities::SP<detran_angle::TabuchiYamamoto>;

//---------------------------------------------------------------------------//
//              end of detran_angle.i
//---------------------------------------------------------------------------//
