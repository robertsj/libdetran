//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_boundary.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran boundary.
 */
//---------------------------------------------------------------------------//

%module(package="detran") boundary
%{
#include <stddef.h>
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryMOC.hh"
#include "boundary/BoundaryCondition.hh"
#include "boundary/FixedBoundary.hh"
%}

// Hide templates from SWIG
%inline
{
#define BOUNDARY_EXPORT
#define BOUNDARY_TEMPLATE_EXPORT(...)
#define BOUNDARY_INSTANTIATE_EXPORT(...)
}

%import "detran_utilities.i"
%import "detran_geometry.i"
%import "detran_angle.i"
%include "transport/DimensionTraits.hh"

%include "BoundaryTraits.hh"
%template(BT1) detran::BoundaryTraits<detran::_1D>;
%template(BT2) detran::BoundaryTraits<detran::_2D>; 
%template(BT3) detran::BoundaryTraits<detran::_3D>; 
%template(BV1) detran::BoundaryValue<detran::_1D>; 
%template(BV2) detran::BoundaryValue<detran::_2D>;
%template(BV3) detran::BoundaryValue<detran::_3D>; 

%include "BoundaryCondition.hh"
%template(BoundCond1D)      detran::BoundaryCondition<detran::_1D>;
%template(BoundCond2D)      detran::BoundaryCondition<detran::_2D>;
%template(BoundCond3D)      detran::BoundaryCondition<detran::_3D>;
%template(BoundCond1DSP)    detran_utilities::SP<detran::BoundaryCondition<detran::_1D> >;
%template(BoundCond2DSP)    detran_utilities::SP<detran::BoundaryCondition<detran::_2D> >;
%template(BoundCond3DSP)    detran_utilities::SP<detran::BoundaryCondition<detran::_3D> >;

%include "FixedBoundary.hh"
// extend fixed bc
//   key = [o, a, group, i, j]
%extend detran_utilities::SP<detran::FixedBoundary<_1D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2])
      value(t, 0, 0, v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2])
      return value(t, 0, 0)
  %}
}
%extend detran_utilities::SP<detran::FixedBoundary<_2D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2])
      value(t, key[3], 0, v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2])
      return value(t, key[3], 0)
  %}
}
%extend detran_utilities::SP<detran::FixedBoundary<_3D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2])
      value(t, key[3], key[4], v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2])
      return value(t, key[3], key[4])
  %}
}
%template(FixedBoundCond1D)     detran::FixedBoundary<detran::_1D>;
%template(FixedBoundCond2D)     detran::FixedBoundary<detran::_2D>;
%template(FixedBoundCond3D)     detran::FixedBoundary<detran::_3D>;
%template(FixedBoundCond1DSP)   detran_utilities::SP<detran::FixedBoundary<detran::_1D> >;
%template(FixedBoundCond2DSP)   detran_utilities::SP<detran::FixedBoundary<detran::_2D> >;
%template(FixedBoundCond3DSP)   detran_utilities::SP<detran::FixedBoundary<detran::_3D> >;

%include "BoundaryBase.hh"
%template(BoundBASE1D)     detran::BoundaryBase<detran::_1D>; 
%template(BoundBASE2D)     detran::BoundaryBase<detran::_2D>;
%template(BoundBASE3D)     detran::BoundaryBase<detran::_3D>; 
%template(BoundBASE1DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_1D> >;
%template(BoundBASE2DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_2D> >;
%template(BoundBASE3DSP)   detran_utilities::SP<detran::BoundaryBase<detran::_3D> >;

%include "BoundaryDiffusion.hh"
// extend diffusion
//   key = [side, inout, group, i, j]
%extend detran_utilities::SP<detran::BoundaryDiffusion<_1D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2])
      value(t, 0, 0, v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2])
      return value(t, 0, 0)
  %}
}
%extend detran_utilities::SP<detran::BoundaryDiffusion<_2D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2])
      value(t, key[3], 0, v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2])
      return value(t, key[3], 0)
  %}
}
%extend detran_utilities::SP<detran::BoundaryDiffusion<_3D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2])
      value(t, key[3], key[4], v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2])
      return value(t, key[3], key[4])
  %}
}
%template(BoundDIFF1D)       detran::BoundaryDiffusion<detran::_1D>;
%template(BoundDIFF1DSP)     detran_utilities::SP<detran::BoundaryDiffusion<detran::_1D> >;
%template(BoundDIFF2D)       detran::BoundaryDiffusion<detran::_2D>;
%template(BoundDIFF2DSP)     detran_utilities::SP<detran::BoundaryDiffusion<detran::_2D> >;
%template(BoundDIFF3D)       detran::BoundaryDiffusion<detran::_3D>;
%template(BoundDIFF3DSP)     detran_utilities::SP<detran::BoundaryDiffusion<detran::_3D> >;

%include "BoundarySN.hh"
// extend sn
//   key = [side, o, a, group, i, j]
%extend detran_utilities::SP<detran::BoundarySN<_1D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2], key[3])
      value(t, 0, 0, v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2], key[3])
      return value(t, 0, 0)
  %}
}
%extend detran_utilities::SP<detran::BoundarySN<_2D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2], key[3])
      value(t, key[4], 0, v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2], key[3])
      return value(t, key[4], 0)
  %}
}
%extend detran_utilities::SP<detran::BoundarySN<_3D> >
{
  %pythoncode 
  %{
    def __setitem__(self, key, v) :
      t = self(key[0], key[1], key[2], key[3])
      value(t, key[4], key[5], v)
    def __getitem__(self, key) :
      t = self(key[0], key[1], key[2], key[3])
      return value(t, key[4], key[5])
  %}
}
%template(BoundSN1D)       detran::BoundarySN<detran::_1D>;
%template(BoundSN1DSP)     detran_utilities::SP<detran::BoundarySN<detran::_1D> >;
%template(BoundSN2D)       detran::BoundarySN<detran::_2D>;
%template(BoundSN2DSP)     detran_utilities::SP<detran::BoundarySN<detran::_2D> >;
%template(BoundSN3D)       detran::BoundarySN<detran::_3D>;
%template(BoundSN3DSP)     detran_utilities::SP<detran::BoundarySN<detran::_3D> >;

%include "BoundaryMOC.hh"
%template(BoundMOC2D)      detran::BoundaryMOC<detran::_2D>;
%template(BoundMOC2DSP)    detran_utilities::SP<detran::BoundaryMOC<detran::_2D> >;

%inline
{

// All the creators return SP's to base classes.  These cast up.
  
// cast base to diffusion
detran_utilities::SP<detran::BoundaryDiffusion<detran::_1D> > 
cast_diff( detran_utilities::SP<detran::BoundaryBase<detran::_1D> >* b )
{
  return detran_utilities::SP<detran::BoundaryDiffusion<detran::_1D> >( *b );
} 
detran_utilities::SP<detran::BoundaryDiffusion<detran::_2D> > 
cast_diff( detran_utilities::SP<detran::BoundaryBase<detran::_2D> >* b )
{
  return detran_utilities::SP<detran::BoundaryDiffusion<detran::_2D> >( *b );
} 
detran_utilities::SP<detran::BoundaryDiffusion<detran::_3D> > 
cast_diff( detran_utilities::SP<detran::BoundaryBase<detran::_3D> >* b )
{
  return detran_utilities::SP<detran::BoundaryDiffusion<detran::_3D> >( *b );
} 

// cast base to sn
detran_utilities::SP<detran::BoundarySN<detran::_1D> > 
cast_sn( detran_utilities::SP<detran::BoundaryBase<detran::_1D> >* b )
{
  return detran_utilities::SP<detran::BoundarySN<detran::_1D> >( *b );
} 
detran_utilities::SP<detran::BoundarySN<detran::_2D> > 
cast_sn( detran_utilities::SP<detran::BoundaryBase<detran::_2D> >* b )
{
  return detran_utilities::SP<detran::BoundarySN<detran::_2D> >( *b );
} 
detran_utilities::SP<detran::BoundarySN<detran::_3D> > 
cast_sn( detran_utilities::SP<detran::BoundaryBase<detran::_3D> >* b )
{
  return detran_utilities::SP<detran::BoundarySN<detran::_3D> >( *b );
} 

// cast base bc to fixed
detran_utilities::SP<detran::FixedBoundary<detran::_1D> > 
cast_fixed( detran_utilities::SP<detran::BoundaryCondition<detran::_1D> >* b )
{
  return detran_utilities::SP<detran::FixedBoundary<detran::_1D> >( *b );
} 
detran_utilities::SP<detran::FixedBoundary<detran::_2D> > 
cast_fixed( detran_utilities::SP<detran::BoundaryCondition<detran::_2D> >* b )
{
  return detran_utilities::SP<detran::FixedBoundary<detran::_2D> >( *b );
} 
detran_utilities::SP<detran::FixedBoundary<detran::_3D> > 
cast_fixed( detran_utilities::SP<detran::BoundaryCondition<detran::_3D> >* b )
{
  return detran_utilities::SP<detran::FixedBoundary<detran::_3D> >( *b );
} 

// return a boundary value
double value(detran::BoundaryTraits<detran::_1D>::value_type* b, int i, int j)           
{ 
  return detran::BoundaryValue<detran::_1D>::value(*b, i, j); 
}
double value(detran::BoundaryTraits<detran::_2D>::value_type* b, int i, int j)           
{ 
  return detran::BoundaryValue<detran::_2D>::value(*b, i, j); 
}
double value(detran::BoundaryTraits<detran::_3D>::value_type* b, int i, int j)           
{ 
  return detran::BoundaryValue<detran::_3D>::value(*b, i, j); 
}
// set a boundary value
void   value(detran::BoundaryTraits<detran::_1D>::value_type* b, int i, int j, double v) 
{ 
  detran::BoundaryValue<detran::_1D>::value(*b, i, j) = v; 
}
void   value(detran::BoundaryTraits<detran::_2D>::value_type* b, int i, int j, double v) 
{ 
  detran::BoundaryValue<detran::_2D>::value(*b, i, j) = v; 
}
void   value(detran::BoundaryTraits<detran::_3D>::value_type* b, int i, int j, double v) 
{ 
  detran::BoundaryValue<detran::_3D>::value(*b, i, j) = v; 
}

} // end inline

//---------------------------------------------------------------------------//
//              end of detran_boundary.i
//---------------------------------------------------------------------------//

