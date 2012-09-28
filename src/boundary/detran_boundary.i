//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_boundary.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran boundary.
 */
//---------------------------------------------------------------------------//

%include "std_vector.i"
%include "BoundaryTraits.hh"
%template(BTrait1D) detran::BoundaryTraits<detran::_1D>;
%template(BTrait2D) detran::BoundaryTraits<detran::_2D>; 
%template(BTrait3D) detran::BoundaryTraits<detran::_3D>; 

%template(BoundaryValue1D)    detran::BoundaryValue<detran::_1D>; 
%template(BoundaryValue2D)    detran::BoundaryValue<detran::_2D>;
%template(BoundaryValue3D)    detran::BoundaryValue<detran::_3D>; 

%include "BoundaryBase.hh"
%template(BoundaryBase1D)     detran::BoundaryBase<detran::_1D>; 
%template(BoundaryBase2D)     detran::BoundaryBase<detran::_2D>;
%template(BoundaryBase3D)     detran::BoundaryBase<detran::_3D>; 
%template(BoundaryBase1DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_1D> >;
%template(BoundaryBase2DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_2D> >;
%template(BoundaryBase3DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_3D> >;

%include "BoundaryDiffusion.hh"
%template(BoundaryDiffusion1D)       detran::BoundaryDiffusion<detran::_1D>;
%template(BoundaryDiffusion1DSP)     detran_utilities::SP<detran::BoundaryDiffusion<detran::_1D> >;
%template(BoundaryDiffusion2D)       detran::BoundaryDiffusion<detran::_2D>;
%template(BoundaryDiffusion2DSP)     detran_utilities::SP<detran::BoundaryDiffusion<detran::_2D> >;
%template(BoundaryDiffusion3D)       detran::BoundaryDiffusion<detran::_3D>;
%template(BoundaryDiffusion3DSP)     detran_utilities::SP<detran::BoundaryDiffusion<detran::_3D> >;

%ignore detran::BoundaryDiffusion<detran::_1D>::operator();

    
%extend detran::BoundaryDiffusion<detran::_1D>
{
  // overload the () so that we can get and set the single double.  using the
  // std_vector functions, the 2D/3D version work just fine.
  double  __getitem__(int s, int g, int io)           { return (*self)(s, g, io); }
  void    __setitem__(int s, int g, int io, double v) { (*self)(s, g, io) = v;    }
}
%extend detran_utilities::SP<detran::BoundaryDiffusion<detran::_1D> >
{
  // overload the () so that we can get and set the single double.  using the
  // std_vector functions, the 2D/3D version work just fine.
  double  __getitem__(int s, int g, int io)           { return (*(*self))(s, g, io); }
  void    __setitem__(int s, int g, int io, double v) { (*(*self))(s, g, io) = v;    }
}
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





