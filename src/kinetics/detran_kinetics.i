//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_kinetics.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran kinetics.
 */
//---------------------------------------------------------------------------//
   
%module(directors="1", allprotected="1", package="detran") kinetics
%{
#include <stddef.h>
// material
#include "kinetics/KineticsMaterial.hh"
#include "kinetics/TimeDependentMaterial.hh"
#include "kinetics/PyTimeDependentMaterial.hh"
#include "kinetics/LinearMaterial.hh"
#include "kinetics/Precursors.hh"
// source
#include "external_source/ExternalSource.hh"
#include "external_source/ConstantSource.hh"
#include "external_source/DiscreteSource.hh"
#include "external_source/IsotropicSource.hh"
#include "kinetics/TimeDependentExternalSource.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/PulsedExternalSource.hh"
//
#include "kinetics/MultiPhysics.hh"
%}

// Hide templates from SWIG
%inline
{
#define KINETICS_EXPORT
#define KINETICS_TEMPLATE_EXPORT(...)
#define KINETICS_INSTANTIATE_EXPORT(...)
}

%include <pycontainer.swg>
%import "detran_utilities.i"
%import "detran_material.i"
%import "detran_angle.i"
%import "detran_geometry.i"
%import "detran_external_source.i"
%import "detran_transport.i"


%include std_string.i
%include callback.i 

// Set the callback setter for python-based time-dependent materials
setCallbackMethod(1, // this is a *unique* identifier
                  detran::PyTimeDependentMaterial, 
                  set_update_impl, 
                  (void* data), (), (), 1)

//// Following construction of a Python material, set the callback function
//%feature("pythonappend") 
//detran::PyTimeDependentMaterial::
//PyTimeDependentMaterial(size_t,size_t,size_t,std::string) 
//%{
//   # the client *must* define update_impl
//   print "PYTHON setting implementation"
//   self.set_update_impl(self.update_material)
//%}

// Sources
%include "TimeDependentExternalSource.hh"
%include "LinearExternalSource.hh"
%include "PulsedExternalSource.hh"
%template(TDSourceSP) detran_utilities::SP<detran::TimeDependentExternalSource>;

// Materials
%include "KineticsMaterial.hh"
%include "TimeDependentMaterial.hh"
%include "PyTimeDependentMaterial.hh"
%include "LinearMaterial.hh"
%include "Precursors.hh"
%template(KineticsMaterialSP)   detran_utilities::SP<detran::KineticsMaterial>;
%template(TDMaterialSP)         detran_utilities::SP<detran::TimeDependentMaterial>;
%template(PYTDMaterialSP)       detran_utilities::SP<detran::PyTimeDependentMaterial>;
%template(vec_material)         std::vector<detran_utilities::SP<detran::KineticsMaterial> >;
%template(PrecursorsSP)         detran_utilities::SP<detran::Precursors>;

// Multiphysics
%include "MultiPhysics.hh"
%template(MultiPhysicsSP)       detran_utilities::SP<detran::MultiPhysics>;

// Downcasts and Upcasts for generic routines
%inline
{
  
  // Kinetics -> Base
  detran_utilities::SP<detran_material::Material> 
  downcast(detran_utilities::SP<detran::KineticsMaterial>* p)
  {
    return detran_utilities::SP<detran_material::Material>(*p);
  } 
  
  // TD -> Base
  detran_utilities::SP<detran_material::Material> 
  downcast(detran_utilities::SP<detran::TimeDependentMaterial>* p)
  {
    return detran_utilities::SP<detran_material::Material>(*p);
  } 
  
  // TD -> PyTD
  detran_utilities::SP<detran::PyTimeDependentMaterial> 
  upcast(detran_utilities::SP<detran::TimeDependentMaterial>* p)
  {
    return detran_utilities::SP<detran::PyTimeDependentMaterial>(*p);
  } 
  
}

%pythoncode
{
  class UserMaterial(object) :
    ''' User-defined time-dependent material.
        
        Users should inherit from this, making sure to
        perform the base construction as well via
          super(UserClass, self).__init__(nm, ng, np, name)
    '''
    
    def __init__(self, nm, ng, np, name="UserTDMaterial") :
      ''' Constructor
    
          Parameters:
            nm(int) -- number of materials
            ng(int) -- number of groups
            np(int) -- number of precursors
      '''
      
      # Create a material smart pointer.  It is returned as a TDMat.
      mat = PyTimeDependentMaterial.Create(nm, ng, np, name)
      self.material = mat
      
      # Set the material update function (the user *must* implement this).
      # Note, we need to upcast to PYTDMat to do so.
      upcast(self.material).set_update_impl(self.update_material)
          
    def update(self, t, dt, order, synthetic) :
      ''' Update the material for this time, step, and order. If needed,
          add the synthetic component.
      ''' 
      self.material.update(t, dt, order, synthetic)
          
    def update_material(self) :
      ''' User-defined update function
      '''
      raise NotImplementedError   
}

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//
