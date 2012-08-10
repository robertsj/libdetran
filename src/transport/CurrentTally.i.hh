//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CurrentTally.i.hh
 * \brief  CurrentTally inline member definitions
 * \author Jeremy Roberts
 * \date   Aug 9, 2012
 */
//---------------------------------------------------------------------------//

#ifndef CURRENTTALLY_I_HH_
#define CURRENTTALLY_I_HH_


namespace detran
{

template <class D>
inline double
CurrentTally<D>::partial_current(const u_int i,
                                 const u_int j,
                                 const u_int k,
                                 const u_int g,
                                 const u_int axis,
                                 const u_int sense)
{
  // Precondition
  Require(axis < D::dimension);
  Require(i <= d_coarsemesh->get_coarse_mesh()->number_cells_x());
  Require(j <= d_coarsemesh->get_coarse_mesh()->number_cells_y());
  Require(k <= d_coarsemesh->get_coarse_mesh()->number_cells_z());
  Require(g < d_number_groups);
  Require(sense < 2);

  u_int nx = d_coarsemesh->get_coarse_mesh()->number_cells_x();
  u_int ny = d_coarsemesh->get_coarse_mesh()->number_cells_y();
  u_int nz = d_coarsemesh->get_coarse_mesh()->number_cells_z();

  u_int index;
  if (axis == 0)
    index = i + j * (nx + 1) + k * (nx + 1) * ny;
  else if (axis == 1)
    index = i + j * nx + k * nx * (ny + 1);
  else
    index = i + j * nx + k * nx * ny;

  return d_partial_current[axis][g][sense][index];
}

template <class D>
inline void
CurrentTally<D>::tally(const u_int i,
                       const u_int j,
                       const u_int k,
                       const u_int g,
                       const u_int o,
                       const u_int a,
                       const face_flux_type psi,
                       const bool incident)
{

  // Make direction triplet
  const u_int dim[] = {i, j, k};

  // Area
  const int area_index[3][2] = { {1, 2}, {0, 2}, {0, 1} };

  // Loop over directions
  for (int d = 0; d < D::dimension; d++)
  {
    // Do we adjust for an incident boundary?
    int inc = 0;
    if (incident) dim[d] ? inc = 1 : inc = -1;

    // Check if on a coarse edge.
    int coarse_edge =
      d_coarsemesh->coarse_edge_flag(d, dim[d] + d_octant_shift[d][o] + inc);
    if (coarse_edge >= 0)
    {

      // Area
      double area = d_coarsemesh->get_fine_mesh()->
                    width(area_index[d][0], dim[area_index[d][0]]) *
                    d_coarsemesh->get_fine_mesh()->
                    width(area_index[d][1], dim[area_index[d][1]]);

      // Tally.
      d_partial_current[d][g][d_octant_shift[d][o]][coarse_edge] +=
        psi[d] * d_quadrature->cosines(d)[a] * d_quadrature->weight(a) * area;

    }
  }
}

// Things simplify a lot for 1D, e.g. the current has no area.
template <>
inline void
CurrentTally<_1D>::tally(const u_int i,
                         const u_int j,
                         const u_int k,
                         const u_int g,
                         const u_int o,
                         const u_int a,
                         const face_flux_type psi,
                         const bool incident)
{

  // Do we adjust for an incident boundary?
  int inc = 0;
  if (incident) i ? inc = 1 : inc = -1;

  // Get the coarse edge index.
  int coarse_edge =
    d_coarsemesh->coarse_edge_flag(i + d_octant_shift[0][o] + inc, 0);

  // If nonnegative, tally the current.
  if (coarse_edge >= 0)
  {
    // Tally.
    d_partial_current[0][g][d_octant_shift[0][o]][coarse_edge] +=
      psi * d_quadrature->mu(0, a) * d_quadrature->weight(a);
  }

}

} // end namespace detran

#endif // CURRENTTALLY_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file CurrentTally.i.hh
//---------------------------------------------------------------------------//
