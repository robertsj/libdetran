//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_transport.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran transport.
 */
//---------------------------------------------------------------------------//

%module detran_transport
%{
// Detran
#include "Acceleration.hh"
#include "Boundary.hh"
//#include "CMR.hh"
#include "FissionSource.hh"
#include "ExternalSource.hh"
#include "ConstantSource.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"
#include "Traits.hh"
// Utilities
#include "Definitions.hh"
#include "SP.hh"
%}

// Load the standard library interfaces.
%include std_vector.i
%include std_string.i
// Standard vec typemaps
%include std_vec_typemap.i

namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

%include "SP.hh"
%include "Boundary.hh"
//%include "FissionSource.hh"
// External/general sources
%include "ExternalSource.hh"
%include "ConstantSource.hh"
//
//%include "CMR.hh"

namespace detran
{

// Simplified State interface.
class State
{
public:
  State(SP<detran::InputDB>       input,
        SP<detran::Mesh>          mesh,
        SP<detran::Quadrature>    quadrature);  
  static SP<State> Create(SP<detran::InputDB>       input,
                          SP<detran::Mesh>          mesh,
                          SP<detran::Quadrature>    quadrature);
  const std::vector<double>& phi(int g) const;
  std::vector<double>& phi(int g);
  const std::vector<double>& psi(int o, int a, int g) const;
  std::vector<double>& psi(int o, int a, int g);
  double eigenvalue() const;
  void set_eigenvalue(double v);
  int number_groups() const;
};
  
// Simplified FissionSource interface.
class FissionSource
{
public:
  FissionSource(SP<detran::State>     state, 
                SP<detran::Mesh>      mesh, 
                SP<detran::Material>  material);
  static SP<FissionSource> Create(SP<detran::State>     state, 
                                  SP<detran::Mesh>      mesh, 
                                  SP<detran::Material>  material);
  void initialize();
  void update();
  void setup_outer(double scale);
  const std::vector<double>& source(int g);
  const std::vector<double>& density();
};

//class Acceleration
//{
//public:
//  Acceleration(SP<detran::Mesh>       mesh, 
//               SP<detran::Material>   material, 
//               SP<detran::Quadrature> quadrature);
//  ~Acceleration(){}
//  virtual void initialize(int level);
//  virtual void tally(int i, int j, int k, int o, int a, int g, double psi) = 0;
//  void homogenize(SP<detran::State> state);
//  int fine_to_coarse(int ijk, int dim);
//  SP<detran::Mesh> get_mesh();
//  SP<detran::Mesh> get_coarse_mesh();
//  SP<detran::Material> get_material();
//  SP<detran::Material> get_coarse_material();
//};

//class CMR : public Acceleration
//{
//public:
//  CMR(SP<detran::Mesh>       mesh, 
//      SP<detran::Material>   material, 
//      SP<detran::Quadrature> quadrature);
//  void tally(int i, int j, int k, int o, int a, int g, double psi);
//};

} // end namespace detran

%inline %{
typedef detran::State::moments_type moments_type;
%}

// Templates 

%template(StateSP)  detran::SP<detran::State>;

%template(FissionSourceSP)  detran::SP<detran::FissionSource>;

%template(ExternalSourceSP) detran::SP<detran::ExternalSource>;
%template(ConstantSourceSP) detran::SP<detran::ConstantSource>;

%template(Boundary1D)    detran::Boundary<detran::_1D>;
%template(Boundary1DSP)  detran::SP<detran::Boundary<detran::_1D> >;
%template(Boundary2D)    detran::Boundary<detran::_2D>;
%template(Boundary2DSP)  detran::SP<detran::Boundary<detran::_2D> >;
%template(Boundary3D)    detran::Boundary<detran::_3D>;
%template(Boundary3DSP)  detran::SP<detran::Boundary<detran::_3D> >;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//





