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

namespace detran
{

//------------------------------------------------------------------------//
// Getters
//------------------------------------------------------------------------//

/// Return total number of cells.
inline int Mesh::number_cells() const
{
  return d_number_cells;
}

/// Return number of cells in specified dimension.
inline int Mesh::number_cells(int dim) const
{
  Require(dim >= 0);
  //Require(dim < d_dimension);
  if (dim == 0)
    return d_number_cells_x;
  else if (dim == 1)
    return d_number_cells_y;
  else
    return d_number_cells_z;
}

inline int Mesh::number_cells_x() const
{
  return d_number_cells_x;
}

inline int Mesh::number_cells_y() const
{
  return d_number_cells_y;
}

inline int Mesh::number_cells_z() const
{
  return d_number_cells_z;
}

inline double Mesh::width(int dim, int ijk) const
{
  Require(dim >= 0);
  Require(dim <  3);
  if (dim == 0)
    return dx(ijk);
  else if (dim == 1)
    return dy(ijk);
  else
    return dz(ijk);
}

inline double Mesh::dx(int i) const
{
  Require (i >= 0);
  Require (i < d_number_cells_x);
  return d_dx[i];
}

inline double Mesh::dy(int j) const
{
  Require (j >= 0);
  Require (j < d_number_cells_y);
  return d_dy[j];
}

inline double Mesh::dz(int k) const
{
  Require (k >= 0);
  Require (k < d_number_cells_z);
  return d_dz[k];
}

inline const vec_dbl& Mesh::dx() const
{
  return d_dx;
}

inline const vec_dbl& Mesh::dy() const
{
  return d_dy;
}

inline const vec_dbl& Mesh::dz() const
{
  return d_dz;
}

inline double Mesh::volume(int cell) const
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

inline int Mesh::dimension() const
{
  return d_dimension;
}


inline int Mesh::index(int i, int j, int k)
{
  Require(i >= 0);
  Require(i < d_number_cells_x);
  Require(j >= 0);
  Require(j < d_number_cells_y);
  Require(k >= 0);
  Require(k < d_number_cells_z);
  return i + j * d_number_cells_x + k * d_number_cells_x * d_number_cells_y;
}

inline int Mesh::cell_to_i(int cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  int i = cell % d_number_cells_x;
  Ensure(i < d_number_cells_x);
  return i;
}

inline int Mesh::cell_to_j(int cell) const
{
  Require(cell >= 0);
  Require(cell < d_number_cells);
  int j = cell % (d_number_cells_x * d_number_cells_y);
  double tmp = std::floor(double(j)/double(d_number_cells_x));
  j = int(tmp);
  Ensure(j < d_number_cells_y);
  return j;
}

inline int Mesh::cell_to_k(int cell) const
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

/// Unimplemented DBC function.
inline bool Mesh::is_valid() const
{
  return true;
}


} // end namespace detran

#endif // MESH_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file Mesh.i.hh
//---------------------------------------------------------------------------//
