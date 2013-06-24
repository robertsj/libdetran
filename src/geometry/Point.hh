//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Point.hh
 *  @brief Point class definition 
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_POINT_HH_
#define detran_geometry_POINT_HH_

#include "geometry/geometry_export.hh"
#include <cmath>
#include <ostream>
#include <iomanip>

namespace detran_geometry
{

/**
 *  @class Point
 *  @brief Represent a point in three-space
 */
class GEOMETRY_EXPORT Point
{

public:

  /// Constructor
  explicit Point(const double xval = 0.0,
                 const double yval = 0.0,
                 const double zval = 0.0);

  Point operator*(const double scale) const
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

  Point operator+(const double p) const
  {
    return Point(d_x+p, d_y+p, d_z+p);
  }

  Point operator-(const double p) const
  {
    return Point(d_x-p, d_y-p, d_z-p);
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

  //@{
  /// Get value for a dimension
  const double& operator[](const size_t dim) const;
  double& operator[](const size_t dim);
  //@}

private:
 
  /// X coordinate
  double d_x;
  /// Y coordinate
  double d_y;
  /// Z coordinate
  double d_z;

};

/// Scale a point.
GEOMETRY_EXPORT Point operator*(const double scale, const Point &p);

/// Distance between two points.
GEOMETRY_EXPORT double distance(const Point &p1, const Point &p2 = Point());

GEOMETRY_EXPORT bool operator> (const Point &P0, const Point &P1);

GEOMETRY_EXPORT bool operator<= (const Point &P0, const Point &P1);

GEOMETRY_EXPORT bool operator< (const Point &P0, const Point &P1);

GEOMETRY_EXPORT bool operator>= (const Point &P0, const Point &P1);

/// Pretty print of a point
GEOMETRY_EXPORT std::ostream& operator<< (std::ostream &out, const Point &p);

/// Pretty print a point
std::ostream& operator<< (std::ostream &out, const Point &p);

} // end namespace detran_geometry

#endif // detran_geometry_POINT_HH_

//----------------------------------------------------------------------------//
//              end of file Point.hh
//----------------------------------------------------------------------------//
