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

namespace detran
{

inline int Acceleration::fine_to_coarse(int ijk, int dim)
{
  Require(dim < b_mesh->dimension());
  Require(ijk >= 0);
  int coarse = 0;
  if (dim == 0)
  {
    Require(ijk < b_mesh->number_cells_x());
    coarse = b_fine_to_coarse_x[ijk];
  }
  else if (dim == 1)
  {
    Require(ijk < b_mesh->number_cells_y());
    coarse = b_fine_to_coarse_y[ijk];
  }
  else if (dim == 2)
  {
    Require(ijk < b_mesh->number_cells_z());
    coarse = b_fine_to_coarse_z[ijk];
  }
  return coarse;
}

} // end namespace detran

#endif /* ACCELERATION_I_HH_ */
