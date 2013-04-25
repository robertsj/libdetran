//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_transport.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran transport.
 */
//---------------------------------------------------------------------------//

%module(package="detran") transport
%{
#include <stddef.h>
#include "transport/DimensionTraits.hh"
#include "transport/FissionSource.hh"
#include "transport/State.hh"
#include "transport/SweepSource.hh"
#include "transport/Homogenize.hh"
%}

// Hide templates from SWIG
%inline
{
#define TRANSPORT_EXPORT
#define TRANSPORT_TEMPLATE_EXPORT(...)
#define TRANSPORT_INSTANTIATE_EXPORT(...)
}

%import "detran_utilities.i"
%import "detran_material.i"
%import "detran_angle.i"
%import "detran_geometry.i"
%import "detran_external_source.i"
%import "detran_boundary.i"


%include "DimensionTraits.hh"

%include "State.hh"
%include "FissionSource.hh"
%include "ScatterSource.hh"
%template(StateSP)          detran_utilities::SP<detran::State>;
%template(FissionSourceSP)  detran_utilities::SP<detran::FissionSource>;
%template(ScatterSourceSP)  detran_utilities::SP<detran::ScatterSource>;

%include "Homogenize.hh"


// We can optionally compile the sweepers again in the future
// if it turns out we want to script some algorithms.

//%include "SweepSource.hh"
//%template(SweepSource1D)    detran::SweepSource<detran::_1D>;
//%template(SweepSource2D)    detran::SweepSource<detran::_2D>;
//%template(SweepSource3D)    detran::SweepSource<detran::_3D>;
//%template(SweepSource1DSP)  detran_utilities::SP<detran::SweepSource<detran::_1D> >;
//%template(SweepSource2DSP)  detran_utilities::SP<detran::SweepSource<detran::_2D> >;
//%template(SweepSource3DSP)  detran_utilities::SP<detran::SweepSource<detran::_3D> >;
//
//// Sweepers
//%include "Sweeper.hh"
//%template(Sweeper1DBASE)     detran::Sweeper<detran::_1D>;
//%template(Sweeper2DBASE)     detran::Sweeper<detran::_2D>;
//%template(Sweeper3DBASE)     detran::Sweeper<detran::_3D>;
//%include "Sweeper1D.hh"
//%include "Sweeper2D.hh"
//%include "Sweeper3D.hh"
//
//%template(Sweeper1D_SD)     detran::Sweeper1D<detran::Equation_SD_1D>;
//%template(Sweeper2D_SD)     detran::Sweeper2D<detran::Equation_SD_2D>;
//%template(Sweeper1D_DD)     detran::Sweeper1D<detran::Equation_DD_1D>;
//%template(Sweeper2D_DD)     detran::Sweeper2D<detran::Equation_DD_2D>;
//%template(Sweeper3D_DD)     detran::Sweeper3D<detran::Equation_DD_3D>;
//%template(Sweeper2D_SC)     detran::Sweeper2D<detran::Equation_SC_2D>;
//%template(Sweeper1D_DDSP)   detran_utilities::SP<detran::Sweeper1D<detran::Equation_DD_1D> >;
//%template(Sweeper2D_DDSP)   detran_utilities::SP<detran::Sweeper2D<detran::Equation_DD_2D> >;
//%template(Sweeper3D_DDSP)   detran_utilities::SP<detran::Sweeper3D<detran::Equation_DD_3D> >;
//%template(Sweeper2D_SCSP)   detran_utilities::SP<detran::Sweeper2D<detran::Equation_SC_2D> >;
//%template(Sweeper1D_SDSP)   detran_utilities::SP<detran::Sweeper1D<detran::Equation_SD_1D> >;
//%template(Sweeper2D_SDSP)   detran_utilities::SP<detran::Sweeper2D<detran::Equation_SD_2D> >;
//
//%include "Sweeper2DMOC.hh"
//%template(Sweeper2DMOC_SC)    detran::Sweeper2DMOC<detran::Equation_SC_MOC>;
//%template(Sweeper2DMOC_SCSP)  detran_utilities::SP<detran::Sweeper2DMOC<detran::Equation_SC_MOC> >;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//
