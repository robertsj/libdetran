//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_moc.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran MOC system.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "MeshMOC.hh"
%include "Point.hh"
%include "PolarQuadrature.hh"
%include "QuadratureMOC.hh"
%include "Segment.hh"
%include "TabuchiYamamoto.hh"
%include "Track.hh"
%include "TrackDB.hh"
%include "Tracker.hh"
%include "Uniform.hh"

%template(MeshMOCSP)          detran::SP<detran::MeshMOC>;
//
%template(QuadratureMOCSP)    detran::SP<detran::QuadratureMOC>;
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
