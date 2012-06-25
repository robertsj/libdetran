//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Segment.hh
 * \brief  Segment class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef SEGMENT_HH_
#define SEGMENT_HH_

// Utilities
#include "DBC.hh"
#include "SP.hh"

// System
#include <iomanip>
#include <ostream>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Segment
 *  \brief One segment along a track
 *
 *  Each segment traverses one flat source region.  To characterize
 *  the segment, we need its length and an index to the region.
 */
//---------------------------------------------------------------------------//

class Segment
{

public:

  typedef SP<Segment> SP_segment;

  Segment(int r, double l)
    : d_region(r)
    , d_length(l)
  {
    Require(r >= 0);
    Require(l >= 0.0);
  }

  int region() const
  {
    return d_region;
  }

  double length() const
  {
    return d_length;
  }

  /// Scale a segment.
  void scale(double v)
  {
    d_length *= v;
  }

private:

  /// \name Private Data
  /// \{

  /// Flat source region this segment crosses.
  int d_region;

  /// Length of this segment.
  double d_length;

  /// \}

};

inline std::ostream& operator<< (std::ostream &out, Segment s)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "(region = " << s.region() << ", length = " << s.length() << ")";
  out.flags(f);
  return out;
}

} // end namespace detran

#endif /* SEGMENT_HH_ */

//---------------------------------------------------------------------------//
//              end of file Segment.hh
//---------------------------------------------------------------------------//
