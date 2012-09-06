//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Mesh.i.hh
 * \brief  Mesh inline member definitions
 * \author Jeremy Roberts
 * \date   Aug 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MESH_I_HH_
#define MESH_I_HH_

namespace detran_geometry
{

//------------------------------------------------------------------------//
// Getters
//------------------------------------------------------------------------//

/// Return total number of cells.
inline Mesh::size_t Mesh::number_cells() const
{
  return d_number_cells;
}

/// Return number of cells in specified dimension.
inline Mesh::size_t Mesh::number_cells(size_t dim) const
{
  if (dim == 0)
    return d_number_cells_x;
  else if (dim == 1)
    return d_number_cells_y;
  else
    return d_number_cells_z;
}

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

inline double Mesh::volume(size_t cell) const
{
  Require(cell < d_number_cells);
  double v = dx(cell_to_i(cell)) *
             dy(cell_to_j(cell)) *
             dz(cell_to_k(cell));
  Ensure(v > 0.0);
  return v;
}

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

inline Mesh::size_t Mesh::dimension() const
{
  return d_dimension;
}


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

inline Mesh::size_t Mesh::cell_to_i(size_t cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  int i = cell % d_number_cells_x;
  Ensure(i < d_number_cells_x);
  return i;
}

inline Mesh::size_t Mesh::cell_to_j(size_t cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  int j = cell % (d_number_cells_x * d_number_cells_y);
  double tmp = std::floor(double(j)/double(d_number_cells_x));
  j = int(tmp);
  Ensure(j < d_number_cells_y);
  return j;
}

inline Mesh::size_t Mesh::cell_to_k(size_t cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  double tmp = std::floor(double(cell)/double(d_number_cells_x*d_number_cells_y));
  int k = int(tmp);
  Ensure(k < d_number_cells_z);
  return k;
}

/// Return a const reference to the full map (useful for IO)
inline const Mesh::mesh_map_type& Mesh::get_mesh_map() const
{
  return d_mesh_map;
}

} // end namespace detran

#endif // MESH_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file Mesh.i.hh
//---------------------------------------------------------------------------//
