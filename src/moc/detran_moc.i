//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_moc.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran MOC system.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

// Includes in rough order of hierarchy
%include "Point.hh"
//
%include "Segment.hh"
%include "Track.hh"
//
%include "PolarQuadrature.hh"
%include "TabuchiYamamoto.hh"
//
%include "QuadratureMOC.hh"
%include "Collocated.hh"
%include "Uniform.hh"
//
%include "TrackDB.hh"
%include "MeshMOC.hh"
%include "Tracker.hh"




%template(MeshMOCSP)          detran::SP<detran::MeshMOC>;
//
%template(QuadratureMOCSP)    detran::SP<detran::QuadratureMOC>;
%template(CollocatedSP)       detran::SP<detran::Collocated>;
%template(UniformSP)          detran::SP<detran::Uniform>;
//
%template(PolarQuadratureSP)  detran::SP<detran::PolarQuadrature>;
%template(TabuchiYamamotoSP)  detran::SP<detran::TabuchiYamamoto>;
//
%template(TrackDBSP)          detran::SP<detran::TrackDB>;
%template(TrackSP)            detran::SP<detran::Track>;
%template(TrackerSP)          detran::SP<detran::Tracker>;

//---------------------------------------------------------------------------//
//              end of detran_moc.i
//---------------------------------------------------------------------------//
