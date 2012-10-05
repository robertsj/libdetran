//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Track.hh
 * \brief  Track class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef TRACK_HH_
#define TRACK_HH_

#include "utilities/Point.hh"
#include "Segment.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <iomanip>
#include <ostream>
#include <vector>

namespace detran_geometry
{

/*!
 *  \class Track
 *  \brief Represents a track across a domain, consisting of several segments
 */
class Track
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Track>         SP_track;
  typedef Segment::SP_segment                 SP_segment;
  typedef detran_utilities::Point             Point;
  typedef std::vector<Segment>                vec_segment;
  typedef vec_segment::const_iterator         iterator;
  typedef vec_segment::const_reverse_iterator riterator;
  typedef detran_utilities::size_t            size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param r0  Entrance point
   *  \param r1   Exit point
   */
  Track(Point r0, Point r1)
    : d_enter(r0)
    , d_exit(r1)
  {
    double d = distance(r0, r1);
    Point p = r1 - r0;
    d_cos_phi = p.x() / d;
    d_sin_phi = p.y() / d;
  }

  /// SP_constructor
  static SP_track Create(Point r0, Point r1)
  {
    SP_track p(new Track(r0, r1));
    return p;
  }

  /// Return my entrance point
  Point enter() const
  {
    return d_enter;
  }

  /// Return my exit point
  Point exit() const
  {
    return d_exit;
  }

  int number_segments() const
  {
    return d_segments.size();
  }

  void add_segment(Segment s)
  {
    d_segments.push_back(s);
  }

  /// Mutable reference to segment
  Segment& segment(size_t i)
  {
    Require(i < d_segments.size());
    return d_segments[i];
  }

  /// Const reference to segment
  const Segment& segment(size_t i) const
  {
    Require(i < d_segments.size());
    return d_segments[i];
  }

  /// Iterator to the beginning of the track
  iterator begin()
  {
    return d_segments.begin();
  }

  /// Iterator to the end of the track
  riterator rbegin()
  {
    return d_segments.rbegin();
  }

  /// Return the track cosine (with respect to x)
  double cos_phi() const
  {
    return d_cos_phi;
  }

  /// Return the track sine (with respect to x)
  double sin_phi() const
  {
    return d_sin_phi;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  friend std::ostream& operator<< (std::ostream &out, Track &t);

};

/// Output stream for a Track
inline std::ostream& operator<< (std::ostream &out, Track &t)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << " enter = " << t.enter()
      << " exit = "  << t.exit() << std::endl;
  for (int s = 0; s < t.number_segments(); s++)
  {
    out << "          segment = " << s << " " << t.segment(s) << std::endl;
  }
  out.flags(f);
  return out;
}

} // end namespace detran_geometry

#endif /* TRACK_HH_ */

//---------------------------------------------------------------------------//
//              end of file Track.hh
//---------------------------------------------------------------------------//
