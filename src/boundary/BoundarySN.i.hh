//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundarySN.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 24, 2012
 *  @brief  Boundary inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYSN_I_HH_
#define detran_BOUNDARYSN_I_HH_

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
inline const typename BoundarySN<D>::bf_type&
BoundarySN<D>::operator()(const size_t side,
                          const size_t o,
                          const size_t a,
                          const size_t g) const
{
  Require(side < D::dimension*2);
  Require(d_quadrature->valid_index(o, a));
  Require(g < d_number_groups);
  size_t angle = d_quadrature->index(o, a);
  return d_boundary_flux[side][g][angle];
}

//---------------------------------------------------------------------------//
template <class D>
inline typename BoundarySN<D>::bf_type&
BoundarySN<D>::operator()(const size_t side,
                          const size_t o,
                          const size_t a,
                          const size_t g)
{
  return const_cast<typename BoundarySN<D>::bf_type&>
  (
    static_cast<const BoundarySN<D>&>(*this)(side, o, a, g)
  );
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundarySN<D>::set(const size_t g)
{
  for(size_t side = 0; side < 2*D::dimension; ++side)
    d_bc[side]->set(g);
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundarySN<D>::update(const size_t g)
{
  for(size_t side = 0; side < 2*D::dimension; ++side)
    d_bc[side]->update(g);
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundarySN<D>::update(const size_t g,
                                  const size_t o,
                                  const size_t a)
{
  for(size_t side = 0; side < 2*D::dimension; ++side)
    d_bc[side]->update(g, o, a);
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundarySN<D>::clear(const size_t g)
{
  for(size_t side = 0; side < 2*D::dimension; ++side)
    for (size_t angle = 0; angle < d_boundary_flux[side][g].size(); ++angle)
      for (size_t j = 0; j < d_boundary_flux[side][g][angle].size(); ++j)
        for (size_t i = 0; i < d_boundary_flux[side][g][angle][0].size(); ++i)
          d_boundary_flux[side][g][angle][j][i] = 0.0;
}

//---------------------------------------------------------------------------//
template <>
inline void BoundarySN<_2D>::clear(const size_t g)
{
  for(size_t side = 0; side < 4; ++side)
    for (size_t angle = 0; angle < d_boundary_flux[side][g].size(); ++angle)
      for (size_t cell = 0; cell < d_boundary_flux[side][g][angle].size(); ++cell)
        d_boundary_flux[side][g][angle][cell] = 0.0;
}

//---------------------------------------------------------------------------//
template <>
inline void BoundarySN<_1D>::clear(const size_t g)
{
  for(size_t side = 0; side < 2; ++side)
    for (size_t angle = 0; angle < d_boundary_flux[side][g].size(); ++angle)
      d_boundary_flux[side][g][angle] = 0.0;
}

//---------------------------------------------------------------------------//
template <class D>
inline void BoundarySN<D>::psi(const size_t  g,
                               double       *v,
                               const int     inout,
                               const int     gs,
                               bool          onlyref)
{
  THROW("NOT IMPLEMENTED");
}

//---------------------------------------------------------------------------//
template <>
inline void BoundarySN<_3D>::psi(const size_t  g,
                                 double       *v,
                                 const int     inout,
                                 const int     gs,
                                 bool          onlyref)
{

  // Loop over all reflective sides and set the flux.
  for (size_t side = 0; side < 6; ++side)
  {
    if ( (!onlyref) || (onlyref && d_is_reflective[side]) )
    {
      for (size_t o = 0; o < 4; ++o)
      {
        size_t octant = d_quadrature->incident_octant(side)[o];
        if (inout == OUT) octant = d_quadrature->outgoing_octant(side)[o];

        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {
          size_t angle = d_quadrature->index(octant, a);
          bf_type &bf = d_boundary_flux[side][g][angle];
          size_t n1 = bf.size();
          size_t n2 = bf[0].size();
          for (size_t i = 0; i < n1; ++i)
          {
            for (size_t j = 0; j < n2; ++j)
            {
              if (gs == SET)
                bf[i][j] = *v;
              else
                *v = bf[i][j];
              v++;
            }
          }
        }
      }
    }
  }
}

//---------------------------------------------------------------------------//
template <>
inline void BoundarySN<_2D>::psi(const size_t  g,
                                 double       *v,
                                 const int     inout,
                                 const int     gs,
                                 bool          onlyref)
{
  // Loop over all reflective sides and set the flux.
  for (size_t side = 0; side < 4; ++side)
  {
    if ( (!onlyref) || (onlyref && d_is_reflective[side]) )
    {
      for (size_t o = 0; o < 2; ++o)
      {
        size_t octant = d_quadrature->incident_octant(side)[o];
        if (inout == OUT) octant = d_quadrature->outgoing_octant(side)[o];

        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {
          size_t angle = d_quadrature->index(octant, a);
          bf_type &bf = d_boundary_flux[side][g][angle];
          size_t n1 = bf.size();
          for (size_t i = 0; i < n1; ++i)
          {
            if (gs == SET)
              bf[i] = *v;
            else
              *v = bf[i];
            v++;
          }
        }
      }
    }
  }
}

//---------------------------------------------------------------------------//
template <>
inline void BoundarySN<_1D>::psi(const size_t  g,
                                 double       *v,
                                 const int     inout,
                                 const int     gs,
                                 bool          onlyref)
{

  // Loop over all reflective sides and set the flux.
  for (size_t side = 0; side < 2; ++side)
  {
    if ( (! onlyref) || (onlyref && d_is_reflective[side]) )
    {
      for (size_t o = 0; o < 1; ++o)
      {
        size_t octant = d_quadrature->incident_octant(side)[o];
        if (inout == OUT) octant = d_quadrature->outgoing_octant(side)[o];

        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {

          size_t angle = d_quadrature->index(octant, a);
          if (gs == SET)
            d_boundary_flux[side][g][angle] = *v;
          else
            *v = d_boundary_flux[side][g][angle];
          v++;
        }
      }
    }
  }
}

} // end namespace detran

#endif /* detran_BOUNDARYSN_I_HH_ */
