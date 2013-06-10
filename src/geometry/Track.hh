//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Track.hh
 *  @brief Track class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_TRACK_HH_
#define detran_geometry_TRACK_HH_

#include "geometry/Point.hh"
#include "geometry/Segment.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <iomanip>
#include <ostream>
#include <vector>

namespace detran_geometry
{

/**
 *  @class Track
 *  @brief Represents a track across a domain, consisting of several segments
 */
/**
 *  @example geometry/test/test_Track.cc
 *
 *  Test of Track class
 */

class GEOMETRY_EXPORT Track
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Track>         SP_track;
  typedef std::vector<Segment>                vec_segment;
  typedef vec_segment::const_iterator         iterator;
  typedef vec_segment::const_reverse_iterator riterator;
  typedef detran_utilities::size_t            size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param r0   entrance point
   *  @param r1   exit point
   *  @param w    track width
   */
  Track(const Point  &r0,
        const Point  &r1,
        const double  w);

  /// SP_constructor
  static SP_track Create(const Point  &r0,
                         const Point  &r1,
                         const double  w);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Add a segment
  void add_segment(const Segment &s);
  /// Mutable reference to segment
  Segment& segment(const size_t i);
  /// Const reference to segment
  const Segment& segment(const size_t i) const;
  /// Return my entrance point
  Point enter() const;
  /// Return my exit point
  Point exit() const;
  /// Iterator to the beginning of the track
  iterator begin();
  /// Iterator to the end of the track
  riterator rbegin();
  /// Return the track cosine
  double cos_phi() const;
  /// Return the track sine
  double sin_phi() const;
  /// Number of segments along this track
  size_t number_segments() const;
  /// Indentifier
  size_t identifier() const;
  /// Return the track width
  double width() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Track segments
  vec_segment d_segments;
  /// Track entrance point
  Point d_enter;
  /// Track exit point
  Point d_exit;
  /// Cosine of the track angle with respect to x
  double d_cos_phi;
  /// Sine of the track angle with respect to y
  double d_sin_phi;
  /// Track width
  double d_width;
  /// Unique identifier
  size_t d_identifier;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  friend std::ostream& operator<< (std::ostream &out, Track &t);

};

/// Output stream for a Track
std::ostream& operator<< (std::ostream &out, Track &t);

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<Track>)
GEOMETRY_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<Track> >)
GEOMETRY_TEMPLATE_EXPORT(std::vector<std::vector<detran_utilities::SP<Track> > >)

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "Track.i.hh"

#endif /* detran_geometry_TRACK_HH_ */

//----------------------------------------------------------------------------//
//              end of file Track.hh
//----------------------------------------------------------------------------//
