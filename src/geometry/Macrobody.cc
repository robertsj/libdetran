//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Macrobody.cc
 *  @brief Macrobody class definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Macrobody.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Macrobody::Macrobody(c_Point &origin)
  : d_origin(origin)
{
  /* ... */
}

//----------------------------------------------------------------------------//
bool Macrobody::contains(c_Point &r) const
{
  bool v = false;
  for (size_t s = 0; s < d_surfaces.size(); ++s)
    v = !d_surfaces[s]->sense(r);
  return v;
}

//----------------------------------------------------------------------------//
Macrobody::vec_point
Macrobody::intersections(const Ray &r, c_dbl t_max)
{
  vec_point body_points, surface_points;
  for (size_t s = 0; s < d_surfaces.size(); ++s)
  {
    surface_points = d_surfaces[s]->intersections(r, t_max);
  }
  return surface_points;
}

//----------------------------------------------------------------------------//
// RPP
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
RightParallelpiped::RightParallelpiped(c_dbl    W,
                                       c_dbl    L,
                                       c_dbl    D,
                                       c_Point &origin)
  : Macrobody(origin)
  , d_width(W)
  , d_length(L)
  , d_depth(D)
{
  /* ... */
}

//----------------------------------------------------------------------------//
bool RightParallelpiped::contains(c_Point &r) const
{
  bool v = false;
  for (size_t s = 0; s < d_surfaces.size(); ++s)
    v = !d_surfaces[s]->sense(r);
  return v;
}

//----------------------------------------------------------------------------//
RightParallelpiped::vec_point
RightParallelpiped::intersections(const Ray &r, c_dbl t_max)
{
  // \todo Implement this
}

//----------------------------------------------------------------------------//
// RCC
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
RightCircularCylinder::
RightCircularCylinder(c_dbl R, c_dbl H, c_Point &origin)
  : Macrobody(origin)
  , d_radius(R)
  , d_height(H)
{
  // \todo Implement this
}

//----------------------------------------------------------------------------//
RightCircularCylinder::vec_point
RightCircularCylinder::intersections(const Ray &r, c_dbl t_max)
{
  // \todo Implement this
}

//----------------------------------------------------------------------------//
// RHP
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
RightHexagonalPrism::
RightHexagonalPrism(c_dbl L, c_dbl H, bool rotate, c_Point &origin)
  : Macrobody(origin)
  , d_length(L)
  , d_height(H)
{
  // \todo Implement this
}

//----------------------------------------------------------------------------//
RightHexagonalPrism::vec_point
RightHexagonalPrism::intersections(const Ray &r, c_dbl t_max)
{
  // \todo Implement this
}

} // namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Macrobody.cc
//----------------------------------------------------------------------------//
