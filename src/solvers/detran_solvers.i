//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_solvers.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%module detran_solvers
%{
// Detran
#include "Boundary.hh"
#include "ExternalSource.hh"
#include "FissionSource.hh"
#include "GaussSeidel.hh"
#include "InnerIteration.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "PowerIteration.hh"
#include "Quadrature.hh"
#include "SourceIteration.hh"
#include "State.hh"
// Utilities
#include "Definitions.hh"
#include "InputDB.hh"
#include "SP.hh"
%}

// Load the standard library interfaces.
%include std_vector.i
%include std_string.i


%include "SP.hh"
%include "InnerIteration.hh"
%include "PowerIteration.hh"

namespace detran
{

// SourceIteration simple SWIG interface.
template <class D>
class SourceIteration: public InnerIteration<D>
{
public:
  // SP Constructor
  static detran::SP<SourceIteration<D> >
  Create(detran::SP<detran::InputDB>         input,
         detran::SP<detran::State>           state,
         detran::SP<detran::Mesh>            mesh,
         detran::SP<detran::Material>        material,
         detran::SP<detran::Quadrature>      quadrature,
         detran::SP<detran::Boundary<D> >    boundary,
         detran::SP<detran::ExternalSource>  q_e,
         detran::SP<detran::FissionSource>   q_f);
  // Solve
  void solve(int g); 
};

// GaussSeidelsimple SWIG interface.
template <class D>
class GaussSeidel
{
public:
  // SP Constructor
  GaussSeidel(detran::SP<detran::InputDB>         input,
              detran::SP<detran::State>           state,
              detran::SP<detran::Mesh>            mesh,
              detran::SP<detran::Material>        material,
              detran::SP<detran::Quadrature>      quadrature,
              detran::SP<detran::Boundary<D> >    boundary,
              detran::SP<detran::ExternalSource>  q_e,
              detran::SP<detran::FissionSource>   q_f);
  // SP Constructor
  static detran::SP<GaussSeidel<D> >
  Create(detran::SP<detran::InputDB>         input,
         detran::SP<detran::State>           state,
         detran::SP<detran::Mesh>            mesh,
         detran::SP<detran::Material>        material,
         detran::SP<detran::Quadrature>      quadrature,
         detran::SP<detran::Boundary<D> >    boundary,
         detran::SP<detran::ExternalSource>  q_e,
         detran::SP<detran::FissionSource>   q_f);
  // Solve
  void solve(); 
};

} // end namespace detran


//%template(StateSP)  detran::SP<detran::State>;

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





