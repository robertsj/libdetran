//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Geometry.cc
 *  @brief Geometry member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "geometry/Geometry.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Geometry::Geometry(const double x, const double y, const double z)
  : d_x(x), d_y(y), d_z(z)
{
  Require(d_x > 0);
  Require(d_y > 0);
  Require(d_z >= 0);
}

//----------------------------------------------------------------------------//
Geometry::SP_geometry
Geometry::Create(const double x, const double y, const double z)
{
  SP_geometry p(new Geometry(x, y, z));
  return p;
}

//----------------------------------------------------------------------------//
void Geometry::add_region(SP_region r)
{
  Require(r);
  d_regions.push_back(r);
  //d_material_map.push_back(m);
}

//----------------------------------------------------------------------------//
Geometry::SP_region Geometry::region(const size_t r)
{
  Require(r < d_regions.size());
  return d_regions[r];
}

//----------------------------------------------------------------------------//
int Geometry::find(const Point &r)
{
  int region = -1;
  for (size_t i = 0; i < number_regions(); i++)
  {
    if (d_regions[i]->contains(r))
    {
      region = i;
      break;
    }
  }
  return region;
}

//----------------------------------------------------------------------------//
Geometry::size_t Geometry::number_regions() const
{
  return d_regions.size();
}

//----------------------------------------------------------------------------//
Geometry::size_t Geometry::material_index(const size_t r) const
{
  Require(r < number_regions());
  return d_material_map[r];
}

} // namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Geometry.cc
//----------------------------------------------------------------------------//

