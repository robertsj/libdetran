//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryMOC.i.hh
 *  @brief  BoundaryMOC inline member definitions.
 *  @author Jeremy Roberts
 *  @date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYMOC_I_HH_
#define detran_BOUNDARYMOC_I_HH_

#include <iostream>

namespace detran
{


//---------------------------------------------------------------------------//
// INHERITED INTERFACE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <class D>
inline void BoundaryMOC<D>::set(const size_t g)
{
  for(int side = 0; side < 2*D::dimension; side++)
    d_bc[side]->set(g);
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundaryMOC<D>::update(const size_t g)
{
  for(int side = 0; side < 2*D::dimension; side++)
    d_bc[side]->update(g);
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundaryMOC<D>::update(const size_t g, const size_t o, const size_t a)
{
  for(int side = 0; side < 2*D::dimension; side++)
    d_bc[side]->update(g, o, a);
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundaryMOC<D>::clear(const size_t g)
{
  /* ... */
}

//---------------------------------------------------------------------------//
// BOUNDARY FLUX ACCESS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <class D>
inline const double&
BoundaryMOC<D>::operator()(const size_t g, const size_t o, const size_t a,
                           const size_t inout, const size_t t) const
{
  Require(d_quadrature->valid_index(o, a));
  Require(g < d_number_groups);
  Require(inout < 2);
  size_t angle = d_quadrature->index(o, a);
  //std::cout << " t = " << t << " size = "
  //          << d_boundary_flux[g][angle][inout].size() << std::endl;
  Require(t < d_boundary_flux[g][angle][inout].size());
  return d_boundary_flux[g][angle][inout][t];
}

//---------------------------------------------------------------------------//
template <class D>
inline double&
BoundaryMOC<D>::operator()(const size_t g, const size_t o, const size_t a,
                           const size_t inout, const size_t t)
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

//---------------------------------------------------------------------------//
template <class D>
inline void
BoundaryMOC<D>::feed_into(const size_t  o1, const size_t  a1, const size_t  t1,
                                size_t &o2,       size_t &a2,       size_t &t2)
{
//  Require(o1 < d_quadrature->number_octants());
//  Require(a1 < d_quadrature->number_angles_octant());
//  int azimuth = d_quadrature->azimuth(a1);
//  Require(t1 < d_quadrature->number_tracks(azimuth));
//  o2 = d_feed_into[o1][a1][t1][0];
//  a2 = d_feed_into[o1][a1][t1][1];
//  t2 = d_feed_into[o1][a1][t1][2];
}

//---------------------------------------------------------------------------//
template <class D>
inline void
BoundaryMOC<D>::feed_from(const size_t  o1, const size_t  a1, const size_t  t1,
                                size_t &o2,       size_t &a2,       size_t &t2)
{
//  Require(o1 < d_quadrature->number_octants());
//  Require(a1 < d_quadrature->number_angles_octant());
//  int azimuth = d_quadrature->azimuth(a1);
//  Require(t1 < d_quadrature->number_tracks(azimuth));
//  o2 = d_feed_from[o1][a1][t1][0];
//  a2 = d_feed_from[o1][a1][t1][1];
//  t2 = d_feed_from[o1][a1][t1][2];
}

} // end namespace detran

#endif // detran_BOUNDARYMOC_I_HH_

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.i.hh
//---------------------------------------------------------------------------//
