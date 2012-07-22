//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryMOC.i.hh
 * \brief  BoundaryMOC inline member definitions.
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYMOC_I_HH_
#define BOUNDARYMOC_I_HH_

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
// INHERITED INTERFACE
//---------------------------------------------------------------------------//

template <class D>
inline void BoundaryMOC<D>::set(int g)
{
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->set(g);
  }
}

template <class D>
inline void BoundaryMOC<D>::update(int g)
{
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->set(g);
  }
}

template <class D>
inline void BoundaryMOC<D>::update(int g, int o, int a)
{
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->update(g, o, a);
  }
}

template <class D>
inline void BoundaryMOC<D>::clear(int g)
{
  /* ... */
}

//---------------------------------------------------------------------------//
// BOUNDARY FLUX ACCESS
//---------------------------------------------------------------------------//

template <class D>
inline const double&
BoundaryMOC<D>::operator()(u_int g, u_int o, u_int a, u_int inout, u_int t) const
{
  Require(d_quadrature->valid_index(o, a));
  Require(g < d_number_groups);
  Require(inout < 2);
  u_int angle = d_quadrature->index(o, a);
  //std::cout << " t = " << t << " size = "
  //          << d_boundary_flux[g][angle][inout].size() << std::endl;
  Require(t < d_boundary_flux[g][angle][inout].size());
  return d_boundary_flux[g][angle][inout][t];
}

template <class D>
inline double&
BoundaryMOC<D>::operator()(u_int g, u_int o, u_int a, u_int inout, u_int t)
{
  // Cast away return type
  return const_cast<double&>
  (
    // Add const to *this's type and call const version
    static_cast<const BoundaryMOC<D>&>(*this)(g, o, a, inout, t)
  );
}

//---------------------------------------------------------------------------//
// INDEXING
//---------------------------------------------------------------------------//

template <class D>
inline void
BoundaryMOC<D>::feed_into(const u_int o1, const u_int a1, const u_int t1,
                          u_int      &o2, u_int      &a2,      u_int &t2)
{
  Require(o1 < d_quadrature->number_octants());
  Require(a1 < d_quadrature->number_angles_octant());
  int azimuth = d_quadrature->azimuth(a1);
  Require(t1 < d_quadrature->number_tracks(azimuth));
  o2 = d_feed_into[o1][a1][t1][0];
  a2 = d_feed_into[o1][a1][t1][1];
  t2 = d_feed_into[o1][a1][t1][2];
}

template <class D>
inline void
BoundaryMOC<D>::feed_from(const u_int o1, const u_int a1, const u_int t1,
                          u_int      &o2, u_int      &a2,      u_int &t2)
{
  Require(o1 < d_quadrature->number_octants());
  Require(a1 < d_quadrature->number_angles_octant());
  int azimuth = d_quadrature->azimuth(a1);
  Require(t1 < d_quadrature->number_tracks(azimuth));
  o2 = d_feed_from[o1][a1][t1][0];
  a2 = d_feed_from[o1][a1][t1][1];
  t2 = d_feed_from[o1][a1][t1][2];
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

// only 2d for now
template class BoundaryMOC<_2D>;


} // end namespace detran

#endif // BOUNDARYMOC_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.i.hh
//---------------------------------------------------------------------------//
