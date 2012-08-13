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

  Assert(index(i, j, k, axis) < d_partial_current[axis][g][sense].size());

  return d_partial_current[axis][g][sense][index(i, j, k, axis)];
}

template <class D>
inline void
CurrentTally<D>::tally(const u_int i,
                       const u_int j,
                       const u_int k,
                       const u_int g,
                       const u_int o,
                       const u_int a,
                       const face_flux_type psi)
{

  // Make direction triplet
  const u_int dim[] = {i, j, k};

  // Area
  const int area_index[3][2] = { {1, 2}, {0, 2}, {0, 1} };

  // Loop over directions
  for (int d = 0; d < D::dimension; d++)
  {
    // Check if on a coarse edge.
    int coarse_edge =
      d_coarsemesh->coarse_edge_flag(dim[d] + d_octant_shift[d][o], d);

    // this coarse index is for the axis of interest.  However, we
    // still need to index into a 2d or 3d map.  We need the IJK of
    // the coarse cell

    if (coarse_edge >= 0)
    {

      // other indices
      int cidx1 = d_coarsemesh->fine_to_coarse(dim[area_index[d][0]], area_index[d][0]);
      int cidx2 = d_coarsemesh->fine_to_coarse(dim[area_index[d][1]], area_index[d][1]);

      int nx = d_coarsemesh->get_coarse_mesh()->number_cells_x();
      int ny = d_coarsemesh->get_coarse_mesh()->number_cells_y();
      int nz = d_coarsemesh->get_coarse_mesh()->number_cells_z();
      int n [] = {nx, ny, nz};

      // Area
      double area = d_coarsemesh->get_fine_mesh()->
                    width(area_index[d][0], dim[area_index[d][0]]) *
                    d_coarsemesh->get_fine_mesh()->
                    width(area_index[d][1], dim[area_index[d][1]]);
      //area = 1.0;

      int idx = 0;
      if (d == 0)
        idx = coarse_edge + (nx + 1) * cidx1 + (nx + 1)*ny * cidx2;
      else if (d == 1)
        idx = cidx1 + coarse_edge * (nx) + (nx)*(ny+1) * cidx2;
      else
        idx = cidx1 + cidx2*(nx) + (nx*ny) * coarse_edge;

      // Tally.
      d_partial_current[d][g][d_octant_shift[d][o]][idx] +=
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
                         const face_flux_type psi)
{
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

// Tallying for incident boundary.
template <class D>
inline void
CurrentTally<D>::tally(const u_int i,
                       const u_int j,
                       const u_int k,
                       const u_int g,
                       const u_int o,
                       const u_int a,
                       const u_int d,
                       const double psi)
{
  using std::cout;
  using std::endl;
  // Make direction triplet
  const u_int dim[] = {i, j, k};
  Require( (dim[d] == 0) or
           (dim[d] == d_coarsemesh->get_fine_mesh()->number_cells(d) - 1) );

  // Area
  const int area_index[3][2] = { {1, 2}, {0, 2}, {0, 1} };

  // Increment for incident direction
  int inc = 0;
  dim[d] ? inc = 1 : inc = -1;

  // Get coarse edge
  int coarse_edge =
    d_coarsemesh->coarse_edge_flag(dim[d] + d_octant_shift[d][o] + inc, d);
  Assert(coarse_edge >= 0);

  // Area
  double area = d_coarsemesh->get_fine_mesh()->
                width(area_index[d][0], dim[area_index[d][0]]) *
                d_coarsemesh->get_fine_mesh()->
                width(area_index[d][1], dim[area_index[d][1]]);
  //area = 1.0;

  // other indices
  int cidx1 = d_coarsemesh->fine_to_coarse(dim[area_index[d][0]], area_index[d][0]);
  int cidx2 = d_coarsemesh->fine_to_coarse(dim[area_index[d][1]], area_index[d][1]);

  int nx = d_coarsemesh->get_coarse_mesh()->number_cells_x();
  int ny = d_coarsemesh->get_coarse_mesh()->number_cells_y();
  int nz = d_coarsemesh->get_coarse_mesh()->number_cells_z();
  int n [] = {nx, ny, nz};

  int idx = 0;
  if (d == 0)
    idx = coarse_edge + (nx + 1) * cidx1 + (nx + 1)*ny * cidx2;
  else if (d == 1)
    idx = cidx1 + coarse_edge * (nx) + (nx)*(ny+1) * cidx2;
  else
    idx = cidx1 + cidx2*(nx) + (nx*ny) * coarse_edge;

  // Tally.
  //cout << " tally ---> d=" << d << " ijk=" << i << j << k << " g="  << g << " sense=" << d_octant_shift[d][o] << " index=" << idx << endl;
  d_partial_current[d][g][d_octant_shift[d][o]][idx] +=
    psi * d_quadrature->cosines(d)[a] * d_quadrature->weight(a) * area;
}

//template <>
//inline void
//CurrentTally<_1D>::tally(const u_int i,
//                         const u_int j,
//                         const u_int k,
//                         const u_int g,
//                         const u_int o,
//                         const u_int a,
//                         const u_int d,
//                         const double psi)
//{
//  // Make direction triplet
//  Require( (i == 0) or
//           (i == d_coarsemesh->get_fine_mesh()->number_cells_x() - 1) );
//
//  // Increment for incident direction
//  int inc = 0;
//  i ? inc = 1 : inc = -1;
//
//  // Get coarse edge
//  int coarse_edge =
//    d_coarsemesh->coarse_edge_flag(i + d_octant_shift[0][o] + inc, 0);
//  Assert(coarse_edge >= 0);
//
//  // Tally.
//  d_partial_current[0][g][d_octant_shift[0][o]][coarse_edge] +=
//    psi * d_quadrature->cosines(0)[a] * d_quadrature->weight(a);
//}

} // end namespace detran

#endif // CURRENTTALLY_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file CurrentTally.i.hh
//---------------------------------------------------------------------------//
