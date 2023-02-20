//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Geometry.cc
 *  @brief Geometry member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "geometry/Geometry.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "geometry/Point.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Geometry::Geometry(double x, double y, double z)
  : d_x(x), d_y(y), d_z(z)
{
  Require(d_x > 0);
  Require(d_y > 0);
  Require(d_z >= 0);
}

//----------------------------------------------------------------------------//
Geometry::Geometry(Mesh::SP_mesh mesh)
  : d_x(mesh->total_width_x()),
    d_y(mesh->total_width_y()),
    d_z(mesh->total_width_z())
{
  Require(d_x > 0);
  Require(d_y > 0);
  Require(d_z >= 0);

  typedef QuadraticSurfaceFactory QSF;

  Require(mesh);
  Require(mesh->dimension() > 1);

  // create a plane-based geometry out of the mesh
  std::vector<double> edges[3];
  Surface::vec_surface surfaces[3];
  for (size_t d = 0; d < 3; ++d)
  {
    size_t n = mesh->number_cells(d);
    edges[d].resize(n + 1);
    surfaces[d].resize(n + 1);
    size_t c[3] = {0, 0, 0};
    c[d] = 1;
    double edge = 0.0;
    for (size_t i = 0; i <= n; ++i)
    {
      edges[d][i]    = edge;
      surfaces[d][i] = QSF::CreatePlane(c[0], c[1], c[2], edge);
      if (i < n) edge += mesh->width(d, i);
    }
  }

  const auto &mat_map = mesh->mesh_map("MATERIAL");
  for (size_t k = 0; k < mesh->number_cells_z(); ++k)
  {
    for (size_t j = 0; j < mesh->number_cells_y(); ++j)
    {
      for (size_t i = 0; i < mesh->number_cells_x(); ++i)
      {
        size_t m    = mat_map[mesh->index(i, j, k)];
        Point b_min, b_max;
        b_min = Point(edges[0][i  ], edges[1][j  ], edges[2][k  ]) - 1e-5;
        b_max = Point(edges[0][i+1], edges[1][j+1], edges[2][k+1]) + 1e-5;
        SP_region r = Region::Create(m, b_min, b_max);
        r->append(surfaces[0][i  ], true);
        r->append(surfaces[0][i+1], false);
        r->append(surfaces[1][j  ], true);
        r->append(surfaces[1][j+1], false);
        if (mesh->dimension() == 3) r->append(surfaces[2][k  ], true);
        if (mesh->dimension() == 3) r->append(surfaces[2][k+1], false);
        add_region(r);
      }
    }
  }

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
}

//----------------------------------------------------------------------------//
Geometry::SP_region Geometry::region(int r) const
{
  Require(r >= 0);
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
int Geometry::number_regions() const
{
  return d_regions.size();
}

//----------------------------------------------------------------------------//
int Geometry::material_index(int r) const
{
  return region(r)->attribute("MATERIAL");
}

} // namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Geometry.cc
//----------------------------------------------------------------------------//

