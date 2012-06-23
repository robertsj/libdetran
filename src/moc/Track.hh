/*
 * Track.hh
 *
 *  Created on: Jun 22, 2012
 *      Author: robertsj
 */

#ifndef TRACK_HH_
#define TRACK_HH_

// Detran
#include "Point.hh"
#include "Segment.hh"

// Utilities
#include "DBC.hh"
#include "SP.hh"

// System
#include <iomanip>
#include <ostream>
#include <vector>


namespace detran
{

class Track: public Object
{

public:

  /// \name Useful Typedefs
  /// \{

  typedef SP<Track>                           SP_track;
  typedef Segment::SP_segment                 SP_segment;
  typedef std::vector<Segment>                vec_segment;
  typedef vec_segment::const_iterator         iterator;
  typedef vec_segment::const_reverse_iterator riterator;

  /// \}

  /*!
   *  \brief Constructor
   *  \param enter  Entrance point
   *  \param exit   Exit point
   */
  Track(Point enter, Point exit)
    : d_enter(enter)
    , d_exit(exit)
  {
    /* ... */
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

  const Segment& segment(int i)
  {
    Require(i < d_segments.size());
    return d_segments[i];
  }

  iterator begin()
  {
    return d_segments.begin();
  }

  riterator rbegin()
  {
    return d_segments.rbegin();
  }

  bool is_valid() const
  {
    return true;
  }

private:

  /// \name Private Data
  /// \{

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

  /// \}

  friend std::ostream& operator<< (std::ostream &out, Track &t);

};

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

} // end namespace detran

#endif /* TRACK_HH_ */

//---------------------------------------------------------------------------//
//              end of file Track.hh
//---------------------------------------------------------------------------//
