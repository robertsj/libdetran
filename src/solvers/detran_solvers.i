//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_solvers.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "InnerIteration.hh"
//%include "SourceIteration.hh"
%include "GaussSeidel.hh"
%include "PowerIteration.hh"

%template(InnerIteration1D)     detran::InnerIteration<detran::_1D>;
%template(InnerIteration1DSP)   detran::SP<detran::InnerIteration<detran::_1D> >;
%template(InnerIteration2D)     detran::InnerIteration<detran::_2D>;
%template(InnerIteration2DSP)   detran::SP<detran::InnerIteration<detran::_2D> >;


%template(GaussSeidel1D)        detran::GaussSeidel<detran::_1D>;
%template(GaussSeidel1DSP)      detran::SP<detran::GaussSeidel<detran::_1D> >;
%template(GaussSeidel2D)        detran::GaussSeidel<detran::_2D>;
%template(GaussSeidel2DSP)      detran::SP<detran::GaussSeidel<detran::_2D> >;

%template(PowerIteration1D)     detran::PowerIteration<detran::_1D>;
%template(PowerIteration1DSP)   detran::SP<detran::PowerIteration<detran::_1D> >;
%template(PowerIteration2D)     detran::PowerIteration<detran::_2D>;
%template(PowerIteration2DSP)   detran::SP<detran::PowerIteration<detran::_2D> >;

namespace detran
{

// This is the one class I couldn't easily wrap.  I narrow it down to it being
// derived and my subsequent use of the "using Base::something" syntax.  SWIG,
// for whatever reason, tries to make getters and setters for those base 
// variables, and the generated code fails due to typedef issues. 
template <class D>
class SourceIteration: public InnerIteration<D>
{
public:
  // Constructor
  SourceIteration(SP<detran::InputDB>          input,
                  SP<detran::State>            state,
                  SP<detran::Mesh>             mesh,
                  SP<detran::Material>         material,
                  SP<detran::Quadrature>       quadrature,
                  SP<detran::BoundaryBase<D> > boundary,
                  SP<detran::ExternalSource>   q_e,
                  SP<detran::FissionSource>    q_f);
  // SP Constructor
  static SP<SourceIteration<D> >
  Create(SP<detran::InputDB>          input,
         SP<detran::State>            state,
         SP<detran::Mesh>             mesh,
         SP<detran::Material>         material,
         SP<detran::Quadrature>       quadrature,
         SP<detran::BoundaryBase<D> > boundary,
         SP<detran::ExternalSource>   q_e,
         SP<detran::FissionSource>    q_f);
  // Solve the within group equation.
  void solve(int g);
};

} // end namespace detran

%template(SourceIteration1D)    detran::SourceIteration<detran::_1D>;
%template(SourceIteration1DSP)  detran::SP<detran::SourceIteration<detran::_1D> >;
%template(SourceIteration2D)    detran::SourceIteration<detran::_2D>;
%template(SourceIteration2DSP)  detran::SP<detran::SourceIteration<detran::_2D> >;


//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//
