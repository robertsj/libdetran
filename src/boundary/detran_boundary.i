//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_boundary.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran boundary.
 */
//---------------------------------------------------------------------------//

%include "BoundaryBase.hh"
%template(BoundaryBase1D)     detran::BoundaryBase<detran::_1D>; 
%template(BoundaryBase2D)     detran::BoundaryBase<detran::_2D>;
%template(BoundaryBase3D)     detran::BoundaryBase<detran::_3D>; 
%template(BoundaryBase1DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_1D> >;
%template(BoundaryBase2DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_2D> >;
%template(BoundaryBase3DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_3D> >;

%include "BoundarySN.hh"
%template(BoundarySN1D)       detran::BoundarySN<detran::_1D>;
%template(BoundarySN1DSP)     detran_utilities::SP<detran::BoundarySN<detran::_1D> >;
%template(BoundarySN2D)       detran::BoundarySN<detran::_2D>;
%template(BoundarySN2DSP)     detran_utilities::SP<detran::BoundarySN<detran::_2D> >;
%template(BoundarySN3D)       detran::BoundarySN<detran::_3D>;
%template(BoundarySN3DSP)     detran_utilities::SP<detran::BoundarySN<detran::_3D> >;

%include "BoundaryMOC.hh"
%template(BoundaryMOC2D)      detran::BoundaryMOC<detran::_2D>;
%template(BoundaryMOC2DSP)    detran_utilities::SP<detran::BoundaryMOC<detran::_2D> >;

//---------------------------------------------------------------------------//
//              end of detran_boundary.i
//---------------------------------------------------------------------------//





