//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CurrentTally.cc
 *  @brief  CurrentTally member definitions
 *  @author Jeremy Roberts
 *  @date   Aug 9, 2012
 */
//---------------------------------------------------------------------------//

#include "CurrentTally.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
CurrentTally<D>::CurrentTally(SP_coarsemesh coarsemesh,
                              SP_quadrature quadrature,
                              const size_t  number_groups)
  : Base(coarsemesh, quadrature, number_groups)
  , d_partial_current(3, vec3_dbl(number_groups,  vec2_dbl(2)))
{
  SP_mesh mesh = d_coarsemesh->get_coarse_mesh();

  size_t nx = d_coarsemesh->get_coarse_mesh()->number_cells_x();
  size_t ny = d_coarsemesh->get_coarse_mesh()->number_cells_y();
  size_t nz = d_coarsemesh->get_coarse_mesh()->number_cells_z();
  size_t n[] = { (nx+1)*ny*nz, nx*(ny+1)*nz, nx*ny*(nz+1) };

  for (int g = 0; g < d_number_groups; g++)
  {
    for (int s = 0; s < 2; s++)
    {
      for (int d = 0; d < D::dimension; d++)
      {
        d_partial_current[d][g][s].resize(n[d], 0.0);
      }
    }
  }

}

//---------------------------------------------------------------------------//
template <class D>
void CurrentTally<D>::reset(const size_t group)
{
  Require(group < d_number_groups);
  for (int d = 0; d < d_partial_current.size(); d++)
    for (int s = 0; s < 2; s++)
      for (int i = 0; i < d_partial_current[d][group][s].size(); i++)
      {
        d_partial_current[d][group][s][i] = 0.0;
      }
}

//---------------------------------------------------------------------------//
template <class D>
void CurrentTally<D>::display()
{
  for (int d = 0; d < d_partial_current.size(); d++)
    for (int g = 0; g < d_partial_current[d].size(); g++)
      for (int s = 0; s < 2; s++)
        for (int i = 0; i < d_partial_current[d][g][s].size(); i++)
        {
          std::cout << " J(d=" << d << ",g=" << g << ",s=" << s << ",i="
                    << i << ") = " << d_partial_current[d][g][s][i]
                    << std::endl;
        }
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
