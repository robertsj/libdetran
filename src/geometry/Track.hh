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
#include "utilities/Iterators.hh"
#include <memory>
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

  typedef std::vector<Segment>                        vec_segment;
  typedef detran_utilities::Reversible<vec_segment>   iterator;
  typedef detran_utilities::size_t                    size_t;
  typedef std::shared_ptr<Track>                      SP_track;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param r0   entrance point
   *  @param r1   exit point
   *  @param w    spatial weight (width in 2d, area in 3d)
   */
  Track(const Point &r0, const Point  &r1, double w);

  SP_track reverse();
  bool reversed() const {return d_reversed;}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Add a segment
  void add_segment(const Segment &s);
  /// Mutable reference to segment
  const Segment& segment(const int i) const
  {
    if (d_reversed)
      return d_segments[number_segments()-i-1];
    return d_segments[i];
  }

  const vec_segment& segments() const { return d_segments;}
  vec_segment& segments() { return d_segments;}

  /// Return my entrance point
  const Point& enter() const;
  /// Return my exit point
  const Point& exit() const;
  /// Iterator to the beginning of the track
  iterator begin();
  /// Iterator to the end of the track
  iterator end();
  /// Number of segments along this track
  int number_segments() const;
  /// Return the track width
  double spatial_weight() const;
  /// Identifying integer
  const int& identifier() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Track segments
  vec_segment& d_segments;
  ///
  vec_segment d_segments_;
  /// Track entrance point
  Point d_enter;
  /// Track exit point
  Point d_exit;
  /// Track spatial weight
  double d_spatial_weight;
  /// Unique identifier
  int d_identifier;
  /// Counts every track added
  static size_t d_number_tracks;
  /// Reversed?
  bool d_reversed;
  /// My reversed one's identifier
  const int& d_reversed_identifier;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  Track(Track &t0, bool reversed);

  friend std::ostream& operator<< (std::ostream &out, Track &t);

  // Comparison by entrance, mu, and xi.  Should be unique.
  bool operator< (const Track &other) const;

  // Equality by entrance, mu, and xi.  Should be unique.
  bool operator== (const Track &other) const;

};

/// Output stream for a Track
std::ostream& operator<< (std::ostream &out, Track &t);

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<Track>)
GEOMETRY_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<Track> >)
GEOMETRY_TEMPLATE_EXPORT(std::vector<std::vector<detran_utilities::SP<Track> > >)

} // end namespace detran_geometry

#endif /* detran_geometry_TRACK_HH_ */

//----------------------------------------------------------------------------//
//              end of file Track.hh
//----------------------------------------------------------------------------//
