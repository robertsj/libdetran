//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Segment.hh
 *  @brief Segment class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_SEGMENT_HH_
#define detran_geometry_SEGMENT_HH_

#include "geometry/geometry_export.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <iomanip>
#include <ostream>
#include <vector>

namespace detran_geometry
{

/**
 *  @class Segment
 *  @brief One segment along a track
 *
 *  Each segment traverses one flat source region.  To characterize
 *  the segment, we need its length and an index to the region.
 *
 *  @todo For acceleration purposes, we might also need to know what
 *        coarse mesh boundary it intersects, if any.
 */
/**
 *  @example geometry/test/test_Segment.cc
 *
 *  Test of Segment class
 */

class GEOMETRY_EXPORT Segment
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::size_t          size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    r   flat source region index
   *  @param    l   segment length
   */
  Segment(const size_t r, const double l);

  /// Return the flat source region
  size_t region() const
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

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Flat source region this segment crosses.
  size_t d_region;
  /// Length of this segment.
  double d_length;

};

/// Output stream for a segment
std::ostream& operator<< (std::ostream &out, const Segment &s);

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<Segment>)
GEOMETRY_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<Segment> >)
GEOMETRY_TEMPLATE_EXPORT(std::vector<Segment>)

} // end namespace detran_geometry

#endif /* detran_geometry_SEGMENT_HH_ */

//----------------------------------------------------------------------------//
//              end of file Segment.hh
//----------------------------------------------------------------------------//
