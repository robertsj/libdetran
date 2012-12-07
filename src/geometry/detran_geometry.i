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

%import "utilities/detran_utilities.i"
%import "angle/detran_angle.i"

%include "Mesh.hh"
//%include "Mesh1D.hh"
//%include "Mesh2D.hh"
//%include "Mesh3D.hh"
%include "PinCell.hh"
%include "Assembly.hh"
%include "Core.hh"
//
%include "MeshMOC.hh"
//%include "Segment.hh"
//%include "Track.hh"
%include "TrackDB.hh"
%include "Tracker.hh"

namespace detran_geometry
{

// Dummy interface to avoid Boost serialization.
class Mesh1D: public Mesh
{
public:  
  Mesh1D(detran_utilities::vec_int xfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_int mat_map);
  Mesh1D(detran_utilities::vec_dbl xfme, 
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_int xfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_dbl xfme, 
         detran_utilities::vec_int mat_map);
};

// Dummy interface to avoid Boost serialization.
class Mesh2D: public Mesh
{
public:
  Mesh2D(detran_utilities::vec_int xfm, 
         detran_utilities::vec_int yfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_dbl ycme, 
         detran_utilities::vec_int mat_map);
  Mesh2D(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_int xfm, 
         detran_utilities::vec_int yfm, 
         detran_utilities::vec_dbl xcme, 
         detran_utilities::vec_dbl ycme, 
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_int mat_map);
};

// Dummy interface to avoid Boost serialization.
class Mesh3D: public Mesh
{
public:
  Mesh3D(detran_utilities::vec_int xfm,  detran_utilities::vec_int yfm,  detran_utilities::vec_int zfm,
         detran_utilities::vec_dbl xcme, detran_utilities::vec_dbl ycme, detran_utilities::vec_dbl zcme,
         detran_utilities::vec_int mat_map);
  Mesh3D(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_dbl zfme, detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_int xfm,  detran_utilities::vec_int yfm,  detran_utilities::vec_int zfm,
         detran_utilities::vec_dbl xcme, detran_utilities::vec_dbl ycme, detran_utilities::vec_dbl zcme,
         detran_utilities::vec_int mat_map);
  static detran_utilities::SP<Mesh>
  Create(detran_utilities::vec_dbl xfme, detran_utilities::vec_dbl yfme, detran_utilities::vec_dbl zfme, detran_utilities::vec_int mat_map);
};

class Segment
{
public:
  typedef detran_utilities::SP<Segment> SP_segment;
  typedef detran_utilities::size_t      size_t;
  Segment(const size_t r, double l);
  int region() const;
  double length() const;
  void scale(double v);
};

class Track
{
public:
  typedef detran_utilities::SP<Track>         SP_track;
  typedef Segment::SP_segment                 SP_segment;
  typedef detran_utilities::Point             Point;
  typedef std::vector<Segment>                vec_segment;
  typedef vec_segment::const_iterator         iterator;
  typedef vec_segment::const_reverse_iterator riterator;
  typedef detran_utilities::size_t            size_t;
  Track(Point r0, Point r1);
  static SP_track Create(Point r0, Point r1);
  Point enter() const;
  Point exit() const;
  int number_segments() const;
  void add_segment(Segment s);
  Segment& segment(size_t i);
  const Segment& segment(size_t i) const;
  iterator begin();
  riterator rbegin();
  double cos_phi() const;
  double sin_phi() const;
};

} // end namespace detran_geometry

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





