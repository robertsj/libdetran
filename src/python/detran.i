//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran library.
 */
//---------------------------------------------------------------------------//

%module pydetran
%{

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
  
// Callow
#include "callow/callow_config.hh"
#include "callow/utils/Initialization.hh"
#include "callow/utils/Typedefs.hh"
//
#include "callow/vector/Vector.hh"
//
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/matrix/MatrixShell.hh"
//
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/Richardson.hh"
#include "callow/solver/Jacobi.hh"
#include "callow/solver/GaussSeidel.hh"
#include "callow/solver/GMRES.hh"
#ifdef CALLOW_ENABLE_PETSC
#include "callow/solver/PetscSolver.hh"
#endif
//
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/preconditioner/PCJacobi.hh"
#include "callow/preconditioner/PCILU0.hh"
#include "callow/preconditioner/PCShell.hh"

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
  
// External source
#include "ExternalSource.hh"
#include "ConstantSource.hh"
#include "DiscreteSource.hh"
#include "IsotropicSource.hh"
  
//---------------------------------------------------------------------------//
// LEVEL 2 
//---------------------------------------------------------------------------//  
  
// Discretization
#include "DimensionTraits.hh"
#include "Equation_DD_1D.hh"
#include "Equation_DD_2D.hh"
#include "Equation_DD_3D.hh"
#include "Equation_SC_2D.hh"
#include "Equation_SD_1D.hh"
#include "Equation_SD_2D.hh"
#include "Equation_SC_MOC.hh"
  
// Boundary
#include "BoundarySN.hh"
#include "BoundaryMOC.hh"
  
// Transport
#include "FissionSource.hh"
#include "State.hh"
#include "SweepSource.hh"
#include "Sweeper.hh"
#include "Sweeper1D.hh"
#include "Sweeper2D.hh"
#include "Sweeper3D.hh"
#include "Sweeper2DMOC.hh"

//---------------------------------------------------------------------------//
// LEVEL 3 
//---------------------------------------------------------------------------//  
      
// Diffusion
#ifdef DETRAN_ENABLE_PETSC
#include "GainOperator.hh"
#include "LossOperator.hh"
#endif
  
// Solvers
#include "InnerIteration.hh"
#include "GaussSeidelMG.hh"
#include "PowerIteration.hh"
#include "SourceIteration.hh"
#include "DiffusionFixedSourceSolver.hh"
#include "DiffusionEigensolver.hh"
#include "DiffusionLossOperator.hh"
#include "DiffusionGainOperator.hh"

//---------------------------------------------------------------------------//
// LEVEL 4 
//---------------------------------------------------------------------------//  

#include "ReactionRates.hh"
#include "Manager.hh"
#include "PyExecute.hh"
#ifdef DETRAN_ENABLE_SILO
#include "SiloOutput.hh"
#endif
#ifdef DETRAN_ENABLE_HDF5
#include "IO_HDF5.hh"
#endif
  
%} // end module pydetran

//------------------------------------//
// CONFIGURATION

%include "detran_config.hh"

//------------------------------------//
// LEVEL 0

// Utilities
%include "detran_utilities.i"

//------------------------------------//
// LEVEL 1

// Callow
%include "callow.i"
// Angle
%include "detran_angle.i"
// Geometry
%include "detran_geometry.i"
// Material
%include "detran_materials.i"
// External source
%include "detran_external_source.i"

//------------------------------------//
// LEVEL 2

// Transport
%include "detran_discretization.i"
%include "detran_boundary.i"
%include "detran_transport.i"

//------------------------------------//
// LEVEL 3

// Diffusion
%include "detran_diffusion.i"

// Transport
%include "detran_solvers.i"

//------------------------------------//
// LEVEL 4

// IO utilities
%include "detran_ioutils.i"

// Post process
%include "detran_postprocess.i"

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
