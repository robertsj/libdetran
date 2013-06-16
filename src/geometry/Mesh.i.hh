//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Mesh.i.hh
 *  @brief Mesh inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_MESH_I_HH_
#define detran_geometry_MESH_I_HH_

namespace detran_geometry
{

//----------------------------------------------------------------------------//
inline Mesh::size_t Mesh::number_cells() const
{
  return d_number_cells;
}

//----------------------------------------------------------------------------//
inline Mesh::size_t Mesh::number_cells(size_t dim) const
{
  if (dim == 0)
    return d_number_cells_x;
  else if (dim == 1)
    return d_number_cells_y;
  else
    return d_number_cells_z;
}

//----------------------------------------------------------------------------//
inline Mesh::size_t Mesh::number_cells_x() const
{
  return d_number_cells_x;
}
inline Mesh::size_t Mesh::number_cells_y() const
{
  return d_number_cells_y;
}
inline Mesh::size_t Mesh::number_cells_z() const
{
  return d_number_cells_z;
}

//----------------------------------------------------------------------------//

inline double Mesh::width(size_t dim, size_t ijk) const
{
  Require(dim <  3);
  if (dim == 0)
    return dx(ijk);
  else if (dim == 1)
    return dy(ijk);
  else
    return dz(ijk);
}
inline double Mesh::dx(size_t i) const
{
  Require (i < d_number_cells_x);
  return d_dx[i];
}
inline double Mesh::dy(size_t j) const
{
  Require (j < d_number_cells_y);
  return d_dy[j];
}
inline double Mesh::dz(size_t k) const
{
  Require (k < d_number_cells_z);
  return d_dz[k];
}

//----------------------------------------------------------------------------//
inline const Mesh::vec_dbl& Mesh::dx() const
{
  return d_dx;
}
inline const Mesh::vec_dbl& Mesh::dy() const
{
  return d_dy;
}
inline const Mesh::vec_dbl& Mesh::dz() const
{
  return d_dz;
}

//----------------------------------------------------------------------------//
inline double Mesh::volume(size_t cell) const
{
  Require(cell < d_number_cells);
  double v = dx(cell_to_i(cell)) *
             dy(cell_to_j(cell)) *
             dz(cell_to_k(cell));
  Ensure(v > 0.0);
  return v;
}

//----------------------------------------------------------------------------//
inline double Mesh::total_width_x() const
{
  return d_total_width_x;
}
inline double Mesh::total_width_y() const
{
  return d_total_width_y;
}
inline double Mesh::total_width_z() const
{
  return d_total_width_z;
}

//----------------------------------------------------------------------------//
inline Mesh::size_t Mesh::dimension() const
{
  return d_dimension;
}

//----------------------------------------------------------------------------//
inline Mesh::size_t Mesh::index(size_t i, size_t j, size_t k)
{
  Require(i >= 0);
  Require(i < d_number_cells_x);
  Require(j >= 0);
  Require(j < d_number_cells_y);
  Require(k >= 0);
  Require(k < d_number_cells_z);
  return i + j * d_number_cells_x + k * d_number_cells_x * d_number_cells_y;
}

//----------------------------------------------------------------------------//
inline Mesh::size_t Mesh::cell_to_i(size_t cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  size_t i = cell % d_number_cells_x;
  Ensure(i < d_number_cells_x);
  return i;
}
inline Mesh::size_t Mesh::cell_to_j(size_t cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  size_t j = cell % (d_number_cells_x * d_number_cells_y);
  double tmp = std::floor(double(j)/double(d_number_cells_x));
  j = size_t(tmp);
  Ensure(j < d_number_cells_y);
  return j;
}
inline Mesh::size_t Mesh::cell_to_k(size_t cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  double tmp = std::floor(double(cell)/double(d_number_cells_x*d_number_cells_y));
  size_t k = size_t(tmp);
  Ensure(k < d_number_cells_z);
  return k;
}

//----------------------------------------------------------------------------//
inline int Mesh::find_cell(Point p)
{
  // Default is -1, meaning not found
  int cell;
  if ( (p.x() < 0.0) || (p.x() > d_total_width_x) ||
       (p.y() < 0.0) || (p.y() > d_total_width_y) ||
       (p.z() < 0.0) || (p.z() > d_total_width_z) )
  {
    cell = -1;
  }
  else
  {
    cell = 0;
    bool found = false;
    size_t ijk[] = {0, 0, 0};
    double xyz[] = {p.x(), p.y(), p.z()};
    for (size_t d = 0; d < dimension(); ++d)
    {
      double current_edge = 0.0;
      for (size_t i = 0; i < number_cells(d); ++i)
      {
        current_edge += width(d, i);
        if (xyz[d] <= current_edge)
        {
          ijk[d] = i;
          found = true;
          break;
        }
      }
    }
    Ensure(found);
    cell = index(ijk[0], ijk[1], ijk[2]);
  }
  return cell;
}

/// Return a const reference to the full map (useful for IO)
inline const Mesh::mesh_map_type& Mesh::get_mesh_map() const
{
  return d_mesh_map;
}

} // end namespace detran

#endif // MESH_I_HH_ 

//----------------------------------------------------------------------------//
//              end of file Mesh.i.hh
//----------------------------------------------------------------------------//
