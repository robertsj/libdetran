//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran library.
 */
//---------------------------------------------------------------------------//

%module pydetran
%{

// Configuration
//#include "config_detran.h"
  
//---------------------------------------------------------------------------//
// LEVEL 0  
//---------------------------------------------------------------------------//
  
// Utilities
#include "Definitions.hh"
#include "DBC.hh"
#include "InputDB.hh"
#include "MathUtilities.hh"
#include "Point.hh"
#include "SP.hh"

//---------------------------------------------------------------------------//
// LEVEL 1
//---------------------------------------------------------------------------//
  
// Angle
#include "Collocated.hh"
#include "GaussLegendre.hh"
#include "LevelSymmetric.hh"
#include "MomentToDiscrete.hh"
#include "PolarQuadrature.hh"
#include "Quadrature.hh"
#include "QuadratureFactory.hh"
#include "QuadratureMOC.hh"
#include "QuadrupleRange.hh"
#include "SphericalHarmonics.hh"
#include "TabuchiYamamoto.hh"
#include "Uniform.hh"
#include "UniformEqual.hh"


// Geometry
#include "Assembly.hh"
#include "Core.hh"
#include "Mesh.hh"
#include "MeshMOC.hh"
#include "Mesh1D.hh" 
#include "Mesh2D.hh" 
#include "Mesh3D.hh"   
#include "PinCell.hh"
#include "Segment.hh"
#include "Track.hh"
#include "TrackDB.hh"
#include "Tracker.hh"

// Material
#include "Material.hh"
  
//---------------------------------------------------------------------------//
// LEVEL 2 
//---------------------------------------------------------------------------//  
  
// Transport
//#include "Acceleration.hh"
#include "BoundaryBase.hh"
#include "Boundary.hh"
#include "BoundaryMOC.hh"
#include "FissionSource.hh"
#include "ExternalSource.hh"
#include "ConstantSource.hh"
#include "DiscreteSource.hh"
#include "State.hh"
#include "SweepSource.hh"
#include "Traits.hh"
//#include "WithinGroupAcceleration.hh"
#include "Sweeper.hh"
#include "Sweeper1D.hh"
#include "Sweeper2D.hh"
#include "Sweeper3D.hh"
#include "Sweeper2DMOC.hh"
//
#include "Equation_DD_1D.hh"
#include "Equation_DD_2D.hh"
#include "Equation_DD_3D.hh"
#include "Equation_SC_2D.hh"
#include "Equation_SD_1D.hh"
#include "Equation_SD_2D.hh"
#include "Equation_SC_MOC.hh"
//
#include "ReactionRates.hh"

//---------------------------------------------------------------------------//
// LEVEL 3 
//---------------------------------------------------------------------------//  
      
// Diffusion
#include "DiffusionEigensolver.hh"
#include "GainOperator.hh"
#include "LossOperator.hh"
  
// Solvers
#include "InnerIteration.hh"
#include "GaussSeidel.hh"
#include "PowerIteration.hh"
#include "SourceIteration.hh"

//---------------------------------------------------------------------------//
// LEVEL 4 
//---------------------------------------------------------------------------//  

#include "Manager.hh"
#include "PyExecute.hh"
  
%} // end module pydetran

//------------------------------------//
// LEVEL 0

// Utilities
%include "detran_utilities.i"

//------------------------------------//
// LEVEL 1

// Angle
%include "detran_angle.i"

// Geometry
%include "detran_geometry.i"

// Material
%include "detran_materials.i"

//------------------------------------//
// LEVEL 2

// Transport
%include "detran_transport.i"

//------------------------------------//
// LEVEL 3

// Diffusion
%include "detran_diffusion.i"

// Transport
%include "detran_solvers.i"

//------------------------------------//
// LEVEL 4

// Anyhere in C/C++ that we need (argc, argv), 

%typemap(in) (int argc, char *argv[]) 
{
  /* Check if is a list */
  if (PyList_Check($input)) 
  {
    int i;
    $1 = PyList_Size($input);
    $2 = (char **) malloc(($1+1)*sizeof(char *));
    for (i = 0; i < $1; i++) 
    {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $2[i] = PyString_AsString(PyList_GetItem($input,i));
      else 
      {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($2);
        return NULL;
      }
    }
    $2[i] = 0;
  } 
  else 
  {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap(freearg) (int argc, char *argv[]) {
  free((char *) $2);
}
%include "PyExecute.hh"
%include "Manager.hh"
%template(Execute1D) detran::PyExecute<detran::_1D>;
%template(Execute2D) detran::PyExecute<detran::_2D>;
%template(Execute3D) detran::PyExecute<detran::_3D>;
