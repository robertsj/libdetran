//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Point.cc
 *  @brief Point class member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "Point.hh"

namespace detran_geometry
{

//---------------------------------------------------------------------------//
Point::Point(const double xval, const double yval, const double zval)
  : d_x(xval)
  , d_y(yval)
  , d_z(zval)
{ 
	/* ... */ 
}

//---------------------------------------------------------------------------//
Point operator*(double scale, const Point &p)
{
  return Point(p.x()*scale, p.y()*scale, p.z()*scale);
}

//---------------------------------------------------------------------------//
double distance(const Point &p1, const Point &p2)
{
  double dx = p2.x()-p1.x();
  double dy = p2.y()-p1.y();
  double dz = p2.z()-p1.z();
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

//---------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream &out, const Point &p)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
  out.flags(f);
  return out;
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Point.cc
//----------------------------------------------------------------------------//
