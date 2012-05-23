//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_solvers.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%include "InnerIteration.hh"
%include "SourceIteration.hh"
%include "GaussSeidel.hh"
%include "PowerIteration.hh"

%template(InnerIteration1D)     detran::InnerIteration<detran::_1D>;
%template(InnerIteration1DSP)   detran::SP<detran::InnerIteration<detran::_1D> >;
%template(InnerIteration2D)     detran::InnerIteration<detran::_2D>;
%template(InnerIteration2DSP)   detran::SP<detran::InnerIteration<detran::_2D> >;

%template(SourceIteration1D)    detran::SourceIteration<detran::_1D>;
%template(SourceIteration1DSP)  detran::SP<detran::SourceIteration<detran::_1D> >;
%template(SourceIteration2D)    detran::SourceIteration<detran::_2D>;
%template(SourceIteration2DSP)  detran::SP<detran::SourceIteration<detran::_2D> >;

%template(GaussSeidel1D)        detran::GaussSeidel<detran::_1D>;
%template(GaussSeidel1DSP)      detran::SP<detran::GaussSeidel<detran::_1D> >;
%template(GaussSeidel2D)        detran::GaussSeidel<detran::_2D>;
%template(GaussSeidel2DSP)      detran::SP<detran::GaussSeidel<detran::_2D> >;

%template(PowerIteration1D)     detran::PowerIteration<detran::_1D>;
%template(PowerIteration1DSP)   detran::SP<detran::PowerIteration<detran::_1D> >;
%template(PowerIteration2D)     detran::PowerIteration<detran::_2D>;
%template(PowerIteration2DSP)   detran::SP<detran::PowerIteration<detran::_2D> >;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//
