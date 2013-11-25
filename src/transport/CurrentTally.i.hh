//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CurrentTally.i.hh
 *  @brief CurrentTally inline member definitions
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_CURRENTTALLY_I_HH_
#define detran_CURRENTTALLY_I_HH_

namespace detran
{

//----------------------------------------------------------------------------//
// Get the partial current on a coarse mesh edge.
template <class D>
inline double
CurrentTally<D>::partial_current(const size_t i,
                                 const size_t j,
                                 const size_t k,
                                 const size_t g,
                                 const size_t axis,
                                 const size_t sense)
{
  Requirev(axis < D::dimension, "Got axis = " + AsString(axis));
  Require(i <= d_coarsemesh->get_coarse_mesh()->number_cells_x());
  Require(j <= d_coarsemesh->get_coarse_mesh()->number_cells_y());
  Require(k <= d_coarsemesh->get_coarse_mesh()->number_cells_z());
  Require(g < d_number_groups);
  Require(sense < 2);

  Assert(index(i, j, k, axis) < d_partial_current[axis][g][sense].size());

  return d_partial_current[axis][g][sense][index(i, j, k, axis)];
}

//----------------------------------------------------------------------------//
// Tallying for main body
template <class D>
inline void
CurrentTally<D>::tally(const size_t i,
                       const size_t j,
                       const size_t k,
                       const size_t g,
                       const size_t o,
                       const size_t a,
                       const face_flux_type psi)
{

  // Make direction triplet
  const size_t dim[] = {i, j, k};

  // Loop over directions of current flow
  for (int d0 = 0; d0 < D::dimension; d0++)
  {

    // Get the coarse edge index along the direction; if not on
    // a coarse edge, -1 is returned.
    int coarse_edge =
      d_coarsemesh->coarse_edge_flag(dim[d0] + d_octant_shift[d0][o], d0);

    if (coarse_edge >= 0)
    {
      // Define directions perpendicular to current
      int d1 = d_perpendicular_index[d0][0];
      int d2 = d_perpendicular_index[d0][1];

      // Coarse mesh triplet
      int cdim[3];
      cdim[d0] = coarse_edge;
      cdim[d1] = d_coarsemesh->fine_to_coarse(dim[d1], d1);
      cdim[d2] = d_coarsemesh->fine_to_coarse(dim[d2], d2);

      // Cardinal coarse edge index
      int idx = index(cdim[0], cdim[1], cdim[2], d0);

      // Area
      double area = d_coarsemesh->get_fine_mesh()->width(d1, dim[d1]) *
                    d_coarsemesh->get_fine_mesh()->width(d2, dim[d2]);

      // Tally
      d_partial_current[d0][g][d_octant_shift[d0][o]][idx] +=
        psi[d0] * d_quadrature->cosines(d0)[a] * d_quadrature->weight(a) * area;

    }
  }
}

//----------------------------------------------------------------------------//
// Tallying for main body, 1D specialization. Things simplify a lot,
// e.g. the current has no area.
template <>
inline void
CurrentTally<_1D>::tally(const size_t i,
                         const size_t j,
                         const size_t k,
                         const size_t g,
                         const size_t o,
                         const size_t a,
                         const face_flux_type psi)
{
  Require(j == 0);
  Require(k == 0);

  // Get the coarse edge index.
  int coarse_edge =
    d_coarsemesh->coarse_edge_flag(i + d_octant_shift[0][o], 0);

  // If nonnegative, tally the current.
  if (coarse_edge >= 0)
  {
    // Tally.
    d_partial_current[0][g][d_octant_shift[0][o]][coarse_edge] +=
      psi * d_quadrature->mu(0, a) * d_quadrature->weight(a);
  }
}

//----------------------------------------------------------------------------//
// Tallying for incident boundary.
template <class D>
inline void
CurrentTally<D>::tally(const size_t i,
                       const size_t j,
                       const size_t k,
                       const size_t g,
                       const size_t o,
                       const size_t a,
                       const size_t d0,
                       const double psi)
{

  // Make direction triplet
  const size_t dim[] = {i, j, k};
  Requirev( (dim[d0] == 0) ||
           (dim[d0] == d_coarsemesh->get_fine_mesh()->number_cells(d0) - 1),
           "Value is " + AsString(dim[d0]));

  // Increment for incident direction
  int inc = 0;
  dim[d0] ? inc = 1 : inc = -1;

  // Get coarse edge
  int coarse_edge =
    d_coarsemesh->coarse_edge_flag(dim[d0] + d_octant_shift[d0][o] + inc, d0);
  Assert(coarse_edge >= 0);

  // Define directions perpendicular to current
  int d1 = d_perpendicular_index[d0][0];
  int d2 = d_perpendicular_index[d0][1];

  // Coarse mesh triplet
  int cdim[3];
  cdim[d0] = coarse_edge;
  cdim[d1] = d_coarsemesh->fine_to_coarse(dim[d1], d1);
  cdim[d2] = d_coarsemesh->fine_to_coarse(dim[d2], d2);

  // Cardinal coarse edge index
  int idx = index(cdim[0], cdim[1], cdim[2], d0);

  // Area
  double area = d_coarsemesh->get_fine_mesh()->width(d1, dim[d1]) *
                d_coarsemesh->get_fine_mesh()->width(d2, dim[d2]);

  // Tally
  d_partial_current[d0][g][d_octant_shift[d0][o]][idx] +=
    psi * d_quadrature->cosines(d0)[a] * d_quadrature->weight(a) * area;
}

} // end namespace detran

#endif // detran_CURRENTTALLY_I_HH_ 

//----------------------------------------------------------------------------//
//              end of file CurrentTally.i.hh
//----------------------------------------------------------------------------//
