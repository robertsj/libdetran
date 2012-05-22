//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Acceleration.i.hh
 * \author Jeremy Roberts
 * \date   May 17, 2012
 * \brief  Acceleration inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef ACCELERATION_I_HH_
#define ACCELERATION_I_HH_

#include <iostream>

namespace detran
{

template <class D>
inline int Acceleration<D>::fine_to_coarse(int ijk, int dim) const
{
  using std::cout;
  using std::endl;

  Require(dim < D::dimension);
  Require(ijk >= 0);
  Require(ijk < b_mesh->number_cells(dim));
  return b_fine_to_coarse[dim][ijk];
}

template <class D>
inline bool Acceleration<D>::on_coarse_boundary(int i, int j, int k, int o) const
{
  bool val = b_coarse_edge_flag[0][i + b_octant_shift[0][o]] *
             b_coarse_edge_flag[1][j + b_octant_shift[1][o]] *
             b_coarse_edge_flag[2][k + b_octant_shift[2][o]];
  return val;
}

template <>
inline bool Acceleration<_2D>::on_coarse_boundary(int i, int j, int k, int o) const
{
  bool val = b_coarse_edge_flag[0][i + b_octant_shift[0][o]] *
             b_coarse_edge_flag[1][j + b_octant_shift[1][o]];
  return val;
}

template <>
inline bool Acceleration<_1D>::on_coarse_boundary(int i, int j, int k, int o) const
{
  bool val = b_coarse_edge_flag[1][i + b_octant_shift[0][o]];
  return val;
}


} // end namespace detran

#endif /* ACCELERATION_I_HH_ */
