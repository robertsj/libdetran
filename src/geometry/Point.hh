//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Point.hh
 *  @brief Point class definition 
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef POINT_HH_
#define POINT_HH_

#include "utilities/utilities_export.hh"
#include <cmath>
#include <ostream>
#include <iomanip>

namespace detran_utilities
{

/**
 *  @class Point
 *  @brief Represent a point in three-space
 */
class UTILITIES_EXPORT Point
{

public:

  /// Constructor
  Point(const double xval = 0.0, 
        const double yval = 0.0,
        const double zval = 0.0);

  Point operator*(double scale) const
  {
    return Point(d_x*scale, d_y*scale, d_z*scale);
  }

  Point operator*(const Point &p) const
  {
    return Point(d_x*p.x(), d_y*p.y(), d_z*p.z());
  }

  Point operator+(const Point &p) const
  {
    return Point(d_x+p.x(), d_y+p.y(), d_z+p.z());
  }

  Point operator-(const Point &p) const
  {
    return Point(d_x-p.x(), d_y-p.y(), d_z-p.z());
  }

  double x() const
  {
    return d_x;
  }
  double y() const
  {
    return d_y;
  }
  double z() const
  {
    return d_z;
  }

private:
 
  /// X coordinate
  double d_x;
  /// Y coordinate
  double d_y;
  /// Z coordinate
  double d_z;

};

/// Scale a point.
UTILITIES_EXPORT Point operator*(double scale, const Point &p);

/// Distance between two points.
UTILITIES_EXPORT double distance(const Point &p1, const Point &p2 = Point());

/// Pretty print of a point
UTILITIES_EXPORT std::ostream& operator<< (std::ostream &out, const Point &p);

/// Pretty print a point
std::ostream& operator<< (std::ostream &out, const Point &p);

} // end namespace detran_utilities

#endif // POINT_HH_ 

//----------------------------------------------------------------------------//
//              end of file Point.hh
//----------------------------------------------------------------------------//
