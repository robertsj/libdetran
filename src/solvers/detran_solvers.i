//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_solvers.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") solvers
%{
#include <stddef.h>
#include "FixedSourceManager.hh"
#include "EigenvalueManager.hh"
#include "time/TimeStepper.hh"
#include "Manager.hh"
#include "time/LRA.hh"
//
#include "kinetics/PyTimeDependentMaterial.hh"
#include "kinetics/LinearMaterial.hh"
%}

%feature("autodoc", "3");

// Hide templates from SWIG
%inline
{
#define SOLVERS_EXPORT
#define SOLVERS_TEMPLATE_EXPORT(...)
#define SOLVERS_INSTANTIATE_EXPORT(...)
}

%import "detran_boundary.i"
%import "detran_kinetics.i"
%import "detran_transport.i"
%include "transport/DimensionTraits.hh"

%include callback.i 

// Set the callback setter for time stepping manager.
// Note, these are indexed 2, 3, and 4
setCallbackMethod(2, 
                  detran::TimeStepper<detran::_1D>, 
                  set_monitor, 
                  (void* data, detran::TimeStepper<detran::_1D>* ts, int step, double t, double dt, int it, bool converged), 
                  (detran::TimeStepper<detran::_1D>* ts, int step, double t, double dt, int it, bool converged), 
                  (ts, step, t, dt, it, converged), 
                  1)
setCallbackMethod(3, 
                  detran::TimeStepper<detran::_2D>, 
                  set_monitor, 
                  (void* data, detran::TimeStepper<detran::_2D>* ts, int step, double t, double dt, int it, bool converged), 
                  (detran::TimeStepper<detran::_2D>* ts, int step, double t, double dt, int it, bool converged), 
                  (ts, step, t, dt, it, converged), 
                  1)
setCallbackMethod(4, 
                  detran::TimeStepper<detran::_3D>, 
                  set_monitor, 
                  (void* data, detran::TimeStepper<detran::_3D>* ts, int step, double t, double dt, int it, bool converged), 
                  (detran::TimeStepper<detran::_3D>* ts, int step, double t, double dt, int it, bool converged), 
                  (ts, step, t, dt, it, converged), 
                  1)
                  
%include "Manager.hh"
%include "TransportManager.hh"                  

%include "FixedSourceManager.hh"
%template(Fixed1D) detran::FixedSourceManager<detran::_1D>;
%template(Fixed2D) detran::FixedSourceManager<detran::_2D>;
%template(Fixed3D) detran::FixedSourceManager<detran::_3D>;

%include "EigenvalueManager.hh"
%template(Eigen1D) detran::EigenvalueManager<detran::_1D>;
%template(Eigen2D) detran::EigenvalueManager<detran::_2D>;
%template(Eigen3D) detran::EigenvalueManager<detran::_3D>;

%include "time/TimeStepper.hh"
%template(Time1D) detran::TimeStepper<detran::_1D>;
%template(Time2D) detran::TimeStepper<detran::_2D>;
%template(Time3D) detran::TimeStepper<detran::_3D>;

%include "time/LRA.hh"
%template(SPLRA) detran_utilities::SP<detran_user::LRA>;

// Downcasts and Upcasts for generic routines
%inline
{
  
  // TD -> LRA
  detran_utilities::SP<detran_user::LRA> 
  as_lra(detran_utilities::SP<detran::TimeDependentMaterial>* p)
  {
    return detran_utilities::SP<detran_user::LRA>(*p);
  } 
  
  // Set LRA physics.  Temporary hack.
  void set_lra_physics(detran::TimeStepper<detran::_2D>* stepper,
                       detran_utilities::SP<detran::MultiPhysics>* physics,
                       detran_utilities::SP<detran::TimeDependentMaterial>* mat)
  {
    stepper->set_multiphysics(*physics,
                              detran_user::update_T_rhs<detran::_2D>,
                              (void *) (*mat).bp());
  }
  
}

//---------------------------------------------------------------------------//
//              end of detran_solvers.i
//---------------------------------------------------------------------------//
