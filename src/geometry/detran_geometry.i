//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_geometry.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran geometry.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") geometry
%{
#include <stddef.h>
#include "utilities/SP.hh"
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
%}

// Hide templates from SWIG
%inline
{
#define GEOMETRY_EXPORT
#define GEOMETRY_TEMPLATE_EXPORT(...)
#define GEOMETRY_INSTANTIATE_EXPORT(...)
#define BOOST_CLASS_EXPORT_KEY(...)
}

%import "utilities/detran_utilities.i"
%import "angle/detran_angle.i"

%include "Mesh.hh"
%include "Mesh1D.hh"
%include "Mesh2D.hh"
%include "Mesh3D.hh"
%include "PinCell.hh"
%include "Assembly.hh"
%include "Core.hh"
//
%include "MeshMOC.hh"
%include "Segment.hh"
%include "Track.hh"
%include "TrackDB.hh"
%include "Tracker.hh"

%template(MeshSP)     detran_utilities::SP<detran_geometry::Mesh>;
%template(Mesh1DSP)   detran_utilities::SP<detran_geometry::Mesh1D>;
%template(Mesh2DSP)   detran_utilities::SP<detran_geometry::Mesh2D>;
%template(Mesh3DSP)   detran_utilities::SP<detran_geometry::Mesh3D>;
%template(PinCellSP)  detran_utilities::SP<detran_geometry::PinCell>;
%template(AssemblySP) detran_utilities::SP<detran_geometry::Assembly>;
%template(CoreSP)     detran_utilities::SP<detran_geometry::Core>;

%template(SegmentSP)  detran_utilities::SP<detran_geometry::Tracker>;
%template(TrackSP)    detran_utilities::SP<detran_geometry::Track>;
%template(TrackDBSP)  detran_utilities::SP<detran_geometry::TrackDB>;
%template(TrackerSP)  detran_utilities::SP<detran_geometry::Tracker>;

//---------------------------------------------------------------------------//
//              end of detran_geometry.i
//---------------------------------------------------------------------------//

