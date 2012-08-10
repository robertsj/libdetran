//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CurrentTally.cc
 * \brief  CurrentTally member definitions
 * \author Jeremy Roberts
 * \date   Aug 9, 2012
 */
//---------------------------------------------------------------------------//

#include "CurrentTally.hh"

namespace detran
{

// Partial currents are sized for groups and sense, which
// should be fine if unused
template <class D>
CurrentTally<D>::CurrentTally(SP_coarsemesh coarsemesh,
                              SP_quadrature quadrature,
                              const u_int number_groups)
  : d_coarsemesh(coarsemesh)
  , d_quadrature(quadrature)
  , d_number_groups(number_groups)
  , d_partial_current(3, vec3_dbl(number_groups,  vec2_dbl(2)))
{
  // Preconditions
  Require(d_coarsemesh);
  Require(d_quadrature);

  SP_mesh mesh = d_coarsemesh->get_coarse_mesh();

  for (int g = 0; g < d_number_groups; g++)
  {
    for (int s = 0; s < 2; s++)
    {
      for (int d = 0; d < D::dimension; d++)
      {
        d_partial_current[d][g][s].resize(mesh->number_cells(d) + 1, 0.0);
      }
    }
  }

  // Octant shift
  d_octant_shift.resize(3, vec_int(8, 0));
  // mu
  d_octant_shift[0][0] =  1;
  d_octant_shift[0][1] =  0;
  d_octant_shift[0][2] =  0;
  d_octant_shift[0][3] =  1;
  d_octant_shift[0][4] =  1;
  d_octant_shift[0][5] =  0;
  d_octant_shift[0][6] =  0;
  d_octant_shift[0][7] =  1;
  // eta
  d_octant_shift[1][0] =  1;
  d_octant_shift[1][1] =  1;
  d_octant_shift[1][2] =  0;
  d_octant_shift[1][3] =  0;
  d_octant_shift[1][4] =  1;
  d_octant_shift[1][5] =  1;
  d_octant_shift[1][6] =  0;
  d_octant_shift[1][7] =  0;
  // xi
  d_octant_shift[2][0] =  1;
  d_octant_shift[2][1] =  1;
  d_octant_shift[2][2] =  1;
  d_octant_shift[2][3] =  1;
  d_octant_shift[2][4] =  0;
  d_octant_shift[2][5] =  0;
  d_octant_shift[2][6] =  0;
  d_octant_shift[2][7] =  0;

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class CurrentTally<_1D>;
template class CurrentTally<_2D>;
template class CurrentTally<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file CurrentTally.cc
//---------------------------------------------------------------------------//
