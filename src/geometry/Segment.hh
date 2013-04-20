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

#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <iomanip>
#include <ostream>

namespace detran_geometry
{

//---------------------------------------------------------------------------//
/*!
 *  \class Segment
 *  \brief One segment along a track
 *
 *  Each segment traverses one flat source region.  To characterize
 *  the segment, we need its length and an index to the region.
 *
 *  \todo For acceleration purposes, we might also need to know what
 *        coarse mesh boundary it intersects, if any.
 */
//---------------------------------------------------------------------------//

class Segment
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Segment> SP_segment;
  typedef detran_utilities::size_t      size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  Segment(const size_t r, double l)
    : d_region(r)
    , d_length(l)
  {
    Require(l >= 0.0);
  }

  /// Return the flat source region
  int region() const
  {
    return d_region;
  }

  /// Return the segment length
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

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Flat source region this segment crosses.
  size_t d_region;

  /// Length of this segment.
  double d_length;

};

/// Output stream for a segment
inline std::ostream& operator<< (std::ostream &out, Segment s)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "(region = " << s.region() << ", length = " << s.length() << ")";
  out.flags(f);
  return out;
}

} // end namespace detran_geometry

#endif /* SEGMENT_HH_ */

//---------------------------------------------------------------------------//
//              end of file Segment.hh
//---------------------------------------------------------------------------//
