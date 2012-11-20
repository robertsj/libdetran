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
#include "callow/solver/EigenSolverCreator.hh"
//
#include "callow/preconditioner/Preconditioner.hh"

// Angle
#include "angle/detran_angle.hh"

// Geometry
#include "geometry/Assembly.hh"
#include "geometry/Core.hh"
#include "geometry/Mesh.hh"
#include "geometry/MeshMOC.hh"
#include "geometry/Mesh1D.hh" 
#include "geometry/Mesh2D.hh" 
#include "geometry/Mesh3D.hh"   
#include "geometry/PinCell.hh"
#include "geometry/Segment.hh"
#include "geometry/Track.hh"
#include "geometry/TrackDB.hh"
#include "geometry/Tracker.hh"

// Material
#include "material/Material.hh"

// External source
#include "external_source/ExternalSource.hh"
#include "external_source/ConstantSource.hh"
#include "external_source/DiscreteSource.hh"
#include "external_source/IsotropicSource.hh"

//---------------------------------------------------------------------------//
// LEVEL 2 
//---------------------------------------------------------------------------//  

// Discretization
#include "discretization/DimensionTraits.hh"
//#include "Equation_DD_1D.hh"
//#include "Equation_DD_2D.hh"
//#include "Equation_DD_3D.hh"
//#include "Equation_SC_2D.hh"
//#include "Equation_SD_1D.hh"
//#include "Equation_SD_2D.hh"
//#include "Equation_SC_MOC.hh"
    
// Boundary
#include "boundary/BoundaryDiffusion.hh"
#include "boundary/BoundarySN.hh"
#include "boundary/BoundaryMOC.hh"

// Transport
#include "transport/FissionSource.hh"
#include "transport/State.hh"
#include "transport/SweepSource.hh"
//#include "Sweeper.hh"
//#include "Sweeper1D.hh"
//#include "Sweeper2D.hh"
//#include "Sweeper3D.hh"
//#include "Sweeper2DMOC.hh"

//---------------------------------------------------------------------------//
// LEVEL 3 
//---------------------------------------------------------------------------//  

// Solvers
#include "FixedSourceManager.hh"
#include "EigenvalueManager.hh"
  
//---------------------------------------------------------------------------//
// LEVEL 4 
//---------------------------------------------------------------------------//  

#include "ReactionRates.hh"
#include "Manager.hh"
//#include "PyExecute.hh"
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
%include "detran_material.i"
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

// Transport
%include "detran_solvers.i"

//------------------------------------//
// LEVEL 4

// IO utilities
%include "detran_ioutils.i"

// Post process
%include "detran_postprocess.i"

//%include "PyExecute.hh"
%include "Manager.hh"
//%template(Execute1D) detran::PyExecute<detran::_1D>;
//%template(Execute2D) detran::PyExecute<detran::_2D>;
//%template(Execute3D) detran::PyExecute<detran::_3D>;
