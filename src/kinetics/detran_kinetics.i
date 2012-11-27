//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_materials.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran kinetics.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"
%include "detran_material.i"
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
%template(vec_source) std::vector<detran_utilities::SP<detran_external_source::ExternalSource> >;

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
  
    def material(self) :
      ''' Return the material object
      '''
      return self.material
          
    def update(self, t, dt, order) :
      ''' Update the material for this time, step, and order
      ''' 
      self.material.update(t, dt, order)
          
    def update_material(self) :
      ''' User-defined update function
      '''
      raise NotImplementedError   
}

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//
