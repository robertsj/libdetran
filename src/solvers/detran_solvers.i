//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_solvers.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%include "detran_kinetics.i"
%include "detran_utilities.i"
%include callback.i 

// Set the callback setter for time stepping manager.
// Note, these are indexed 2, 3, and 4
setCallbackMethod(2, 
                  detran::TimeStepper<detran::_1D>, 
                  set_monitor, 
                  (void* data, detran::TimeStepper<detran::_1D>* ts, int step, double t, double dt, int it), 
                  (detran::TimeStepper<detran::_1D>* ts, int step, double t, double dt, int it), 
                  (ts, step, t, dt, it), 
                  1)
setCallbackMethod(3, 
                  detran::TimeStepper<detran::_2D>, 
                  set_monitor, 
                  (void* data, detran::TimeStepper<detran::_2D>* ts, int step, double t, double dt, int it), 
                  (detran::TimeStepper<detran::_2D>* ts, int step, double t, double dt, int it), 
                  (ts, step, t, dt, it), 
                  1)
setCallbackMethod(4, 
                  detran::TimeStepper<detran::_3D>, 
                  set_monitor, 
                  (void* data, detran::TimeStepper<detran::_3D>* ts, int step, double t, double dt, int it), 
                  (detran::TimeStepper<detran::_3D>* ts, int step, double t, double dt, int it), 
                  (ts, step, t, dt, it), 
                  1)
                  
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

//---------------------------------------------------------------------------//
//              end of detran_solvers.i
//---------------------------------------------------------------------------//
