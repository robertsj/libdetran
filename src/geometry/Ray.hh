//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Ray.hh
 *  @brief Ray class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_RAY_HH_
#define detran_geometry_RAY_HH_

#include "geometry/Point.hh"
#include "utilities/SoftEquivalence.hh"
#include "utilities/TinyVector.hh"

namespace detran_geometry
{

/// Defines a three-dimensional ray
struct Ray
{
  typedef detran_utilities::TinyVector<unsigned int, 3>   tiny_3;
  Ray(const Point &O, const Point &D)
    : origin(O)
    , direction(D)
    , inv_direction(1.0/direction.x(), 1.0/direction.y(), 1.0/direction.z())
    , sign(inv_direction.x() < 0, inv_direction.y() < 0, inv_direction.z() < 0)
  {
    Require(detran_utilities::soft_equiv(distance(direction), 1.0));
  }

  /// Origin of the ray
  const Point origin;
  /// Direction of the ray
  const Point direction;
  /// Elementwise-inverse of the direction
  const Point inv_direction;
  /// Sign of each direction component
  const tiny_3 sign;

};

inline std::ostream& operator<< (std::ostream &out, const Ray &p)
{
  out << " origin=" << p.origin <<  " direction=" << p.direction;
  return out;
}

} // end namespace detran_geometry

#endif /* detran_geometry_RAY_HH_ */

//----------------------------------------------------------------------------//
//              end of file Ray.hh
//----------------------------------------------------------------------------//
