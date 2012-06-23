//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Point.hh
 * \brief  Point 
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef POINT_HH_
#define POINT_HH_

// System
#include <cmath>
#include <ostream>
#include <iomanip>

namespace detran
{

class Point
{

public:


  Point(double xval = 0.0, double yval = 0.0)
    : d_x(xval)
    , d_y(yval)
  { /* ... */ }

  Point operator*(double scale) const
  {
    return Point(d_x*scale, d_y*scale);
  }

  Point operator*(Point p) const
  {
    return Point(d_x*p.x(), d_y*p.y());
  }

  Point operator+(Point p) const
  {
    return Point(d_x+p.x(), d_y+p.y());
  }

  Point operator-(Point p) const
  {
    return Point(d_x-p.x(), d_y-p.y());
  }

  double x() const
  {
    return d_x;
  }
  double y() const
  {
    return d_y;
  }

private:

  /// \name Private Data
  /// \{

  double d_x;
  double d_y;

  /// \}

};

/// Scale a point.
inline Point operator*(double scale, Point p)
{
  return Point(p.x()*scale, p.y()*scale);
}

/// Distance between two points.
inline double distance(Point p1, Point p2)
{
  double dx = p2.x()-p1.x();
  double dy = p2.y()-p1.y();
  return std::sqrt(dx*dx + dy*dy);
}

inline std::ostream& operator<< (std::ostream &out, Point p)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "(" << p.x() << ", " << p.y() << ")";
  out.flags(f);
  return out;
}

} // end namespace detran

#endif // POINT_HH_ 

//---------------------------------------------------------------------------//
//              end of file Point.hh
//---------------------------------------------------------------------------//
