//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GeometryMesher.cc
 *  @brief GeometryMesher
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "geometry/GeometryMesher.hh"
#include "utilities/MathUtilities.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
GeometryMesher::GeometryMesher(const double max_spacing,
                               const size_t max_cells)
  : d_max_spacing(max_spacing)
  , d_max_cells(max_cells)
{

}

//----------------------------------------------------------------------------//
void GeometryMesher::meshify(SP_geometry geo, const size_t scheme)
{
  if (geo->dimension() == 2)
  {
    meshify_2D(geo, scheme);
  }
}

//----------------------------------------------------------------------------//
void GeometryMesher::meshify_2D(SP_geometry geo, const size_t scheme)
{
  double W_X = geo->width_x();
  if (scheme == CONSTANT)
  {
    size_t nx = std::min(geo->width_x() / d_max_spacing, (double)d_max_cells);
    size_t ny = std::min(geo->width_y() / d_max_spacing, (double)d_max_cells);
    d_x_edges = detran_utilities::linspace(0.0, geo->width_x(), nx);
    d_y_edges = detran_utilities::linspace(0.0, geo->width_y(), ny);
  }
  else
  {

  }
}


} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Geometry.hh
//----------------------------------------------------------------------------//

