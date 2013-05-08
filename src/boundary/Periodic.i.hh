//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Periodic.i.hh
 *  @brief Periodic inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//


#ifndef detran_PERIODIC_I_HH_
#define detran_PERIODIC_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
void Periodic<D>::update(const size_t g)
{
  for (size_t o = 0; o < d_quadrature->number_octants()/2; ++o)
  {
    for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
    {
      (*d_boundary)(d_periodic_side, d_octants[o][Boundary_T::IN], a, g) =
        (*d_boundary)(d_side, d_octants[o][Boundary_T::OUT], a, g);
    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
void Periodic<D>::update(const size_t g, const size_t o, const size_t a)
{
  // Determine if an outgoing octant, and if so, update the periodic surface
  for (size_t i = 0; i < d_octants.size(); ++i)
    if (d_octants[i][Boundary_T::OUT] == o)
      (*d_boundary)(d_periodic_side, o, a, g) = (*d_boundary)(d_side, o, a, g);
}

BOUNDARY_INSTANTIATE_EXPORT(Periodic<_1D>)
BOUNDARY_INSTANTIATE_EXPORT(Periodic<_2D>)
BOUNDARY_INSTANTIATE_EXPORT(Periodic<_3D>)

} // end namespace periodic

#endif /* detran_PERIODIC_I_HH_ */
