//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Point.cc
 *  @brief Point class member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "geometry/Point.hh"
#include "utilities/DBC.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
// POINT MEMBERS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
Point::Point(const double xval, const double yval, const double zval)
  : d_x(xval)
  , d_y(yval)
  , d_z(zval)
{ 
	/* ... */ 
}

//----------------------------------------------------------------------------//
const double& Point::operator[](const size_t dim) const
{
  Require(dim < 3);
  return dim == 0 ? d_x : (dim == 1 ? d_y : d_z);
}

//----------------------------------------------------------------------------//
double& Point::operator[](const size_t dim)
{
  Require(dim < 3);
  return dim == 0 ? d_x : (dim == 1 ? d_y : d_z);
}

//----------------------------------------------------------------------------//
// POINT HELPERS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
Point operator*(const double scale, const Point &p)
{
  return Point(p.x()*scale, p.y()*scale, p.z()*scale);
}

//----------------------------------------------------------------------------//
double distance(const Point &p1, const Point &p2)
{
  double dx = p2.x()-p1.x();
  double dy = p2.y()-p1.y();
  double dz = p2.z()-p1.z();
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

//----------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream &out, const Point &p)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
  out.flags(f);
  return out;
}

//----------------------------------------------------------------------------//
bool operator< (const Point &P0, const Point &P1)
{
  return (P0.x() < P1.x()) && (P0.y() < P1.y()) && (P0.z() < P1.z());
}

//----------------------------------------------------------------------------//
bool operator<= (const Point &P0, const Point &P1)
{
  return (P0.x() <= P1.x()) && (P0.y() <= P1.y()) && (P0.z() <= P1.z());
}

//----------------------------------------------------------------------------//
bool operator> (const Point &P0, const Point &P1)
{
  return (P0.x() > P1.x()) && (P0.y() > P1.y()) && (P0.z() > P1.z());
}

//----------------------------------------------------------------------------//
bool operator>= (const Point &P0, const Point &P1)
{
  return (P0.x() >= P1.x()) && (P0.y() >= P1.y()) && (P0.z() >= P1.z());
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Point.cc
//----------------------------------------------------------------------------//
