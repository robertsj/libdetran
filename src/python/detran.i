//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran library.
 */
//---------------------------------------------------------------------------//
namespace
{
}
%module pydetran
%{

//---------------------------------------------------------------------------//
// LEVEL 0  
//---------------------------------------------------------------------------//
  
// Utilities
#include "Definitions.hh"
#include "DBC.hh"
#include "InputDB.hh"
#include "SP.hh"

//---------------------------------------------------------------------------//
// LEVEL 1
//---------------------------------------------------------------------------//
  
// Angle
#include "GaussLegendre.hh"
#include "LevelSymmetric.hh"
#include "MomentToDiscrete.hh"
#include "Quadrature.hh"
#include "QuadrupleRange.hh"
#include "SphericalHarmonics.hh"
#include "UniformEqual.hh"

// Geometry
#include "Mesh.hh"
#include "Mesh1D.hh" 
#include "Mesh2D.hh" 
#include "Mesh3D.hh"   
#include "PinCell.hh"
#include "Assembly.hh"
#include "Core.hh"
  
// Material
#include "Material.hh"
  
//---------------------------------------------------------------------------//
// LEVEL 2 
//---------------------------------------------------------------------------//  
  
// Transport
//#include "Acceleration.hh"
#include "Boundary.hh"
#include "FissionSource.hh"
#include "ExternalSource.hh"
#include "ConstantSource.hh"
#include "State.hh"
#include "SweepSource.hh"
#include "Traits.hh"
//#include "WithinGroupAcceleration.hh"
#include "Sweeper.hh"
#include "Sweeper1D.hh"
#include "Sweeper2D.hh"
#include "Sweeper3D.hh"
//
#include "Equation_DD_1D.hh"
#include "Equation_DD_2D.hh"
#include "Equation_DD_3D.hh"
#include "Equation_SC_2D.hh"

//---------------------------------------------------------------------------//
// LEVEL 3 
//---------------------------------------------------------------------------//  
      
// Solvers
#include "InnerIteration.hh"
#include "GaussSeidel.hh"
#include "PowerIteration.hh"
#include "SourceIteration.hh"

//---------------------------------------------------------------------------//
// LEVEL 4 
//---------------------------------------------------------------------------//  

// TBD (manager, post process, etc.)
  
%}

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

// Transport
%include "detran_solvers.i"

