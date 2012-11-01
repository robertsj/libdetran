//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundarySN.i.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Boundary inline member definitions.
 */
//---------------------------------------------------------------------------//
#ifndef BOUNDARYSN_I_HH_
#define BOUNDARYSN_I_HH_

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
// Access via cardinal indices (i.e. for sweeping)
//---------------------------------------------------------------------------//

template <class D>
inline const typename BoundarySN<D>::boundary_flux_type&
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

template <class D>
inline typename BoundarySN<D>::boundary_flux_type&
BoundarySN<D>::operator()(const size_t side,
                          const size_t o,
                          const size_t a,
                          const size_t g)
{
  // Cast away return type
  return const_cast<typename BoundarySN<D>::boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const BoundarySN<D>&>(*this)(side, o, a, g)
  );
}

//---------------------------------------------------------------------------//
// Access to incident and outgoing flux with local ordering
//---------------------------------------------------------------------------//

template <class D>
inline const typename BoundarySN<D>::boundary_flux_type&
BoundarySN<D>::incident(const size_t side,
                        const size_t angle,
                        const size_t g) const
{
  Require(side >= 0);
  Require(side < D::dimension*2);
  Require(angle >= 0);
  Require(angle < d_quadrature->number_angles_octant());
  Require(g >= 0);
  Require(g < d_number_groups);
  //return d_boundary_flux[side][g][angle];
  int index = ordered_angle(side, angle, IN);
  return d_boundary_flux[side][g][index];
}

template <class D>
inline typename BoundarySN<D>::boundary_flux_type&
BoundarySN<D>::incident(const size_t side,
                        const size_t angle,
                        const size_t g)
{
  // Cast away return type
  return const_cast<typename BoundarySN<D>::boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const BoundarySN<D>&>(*this).incident(side, angle, g)
  );
}

template <class D>
inline const typename BoundarySN<D>::boundary_flux_type&
BoundarySN<D>::outgoing(const size_t side, const size_t angle, const size_t g) const
{
  int index = ordered_angle(side, angle, OUT);
}

template <class D>
inline typename BoundarySN<D>::boundary_flux_type&
BoundarySN<D>::outgoing(const size_t side, const size_t angle, const size_t g)
{
  // Cast away return type
  return const_cast<boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const BoundarySN<D>&>(*this).outgoing(side, angle, g)
  );
}

template <class D>
inline void BoundarySN<D>::set(const size_t g)
{
  // Set boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->set(g);
  }
}

template <class D>
inline void BoundarySN<D>::update(const size_t g)
{
  // Set boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->update(g);
  }
}

template <class D>
inline void BoundarySN<D>::update(const size_t g,
                                  const size_t o,
                                  const size_t a)
{
  // Set boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->update(g, o, a);
  }
}

template <class D>
inline void BoundarySN<D>::clear(const size_t g)
{
  // Clear the boundary flux.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    for (int angle = 0; angle < d_boundary_flux[side][g].size(); angle++)
    {
      for (int j = 0; j < d_boundary_flux[side][g][angle].size(); j++)
      {
        for (int i = 0; i < d_boundary_flux[side][g][angle][0].size(); i++)
        {
          d_boundary_flux[side][g][angle][j][i] = 0.0;
        }
      }
    }
  }
}

template <>
inline void BoundarySN<_2D>::clear(const size_t g)
{
  // Clear the boundary flux.
  for(int side = 0; side < 4; side++)
  {
    for (int angle = 0; angle < d_boundary_flux[side][g].size(); angle++)
    {
      for (int cell = 0; cell < d_boundary_flux[side][g][angle].size(); cell++)
      {
        d_boundary_flux[side][g][angle][cell] = 0.0;
      }
    }
  }
}

template <>
inline void BoundarySN<_1D>::clear(const size_t g)
{
  // Clear the boundary flux.
  for(int side = 0; side < 2; side++)
  {
    for (int angle = 0; angle < d_boundary_flux[side][g].size(); angle++)
    {
      d_boundary_flux[side][g][angle] = 0.0;
    }
  }
}

//---------------------------------------------------------------------------//
// Local ordering index.
//---------------------------------------------------------------------------//

template <class D>
inline detran_utilities::size_t
BoundarySN<D>::ordered_angle(const size_t side,
                             const size_t angle,
                             const size_t inout) const
{
  return 0;
}

template <>
inline detran_utilities::size_t
BoundarySN<_2D>::ordered_angle(const size_t side,
                               const size_t angle,
                               const size_t inout) const
{
  return 0;
}

template <>
inline detran_utilities::size_t
BoundarySN<_1D>::ordered_angle(const size_t side,
                               const size_t angle,
                               const size_t inout) const
{
  if (inout == IN)
  {

  }
}

// Incident condition setter and getter for Krylov solvers.

template <class D>
inline void BoundarySN<D>::psi(const size_t g,
                               double *v,
                               const int inout,
                               const int gs,
                               bool onlyref)
{
  THROW("NOT IMPLEMENTED");
}

template <>
inline void BoundarySN<_3D>::psi(const size_t g,
                                 double *v,
                                 const int inout,
                                 const int gs,
                                 bool onlyref)
{

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 6; side++)
  {
    if ( (not onlyref) or (onlyref and d_is_reflective[side]) )
    {
      for (int o = 0; o < 4; o++)
      {
        int octant = d_quadrature->incident_octant(side)[o];
        if (inout == OUT) octant = d_quadrature->outgoing_octant(side)[o];

        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(octant, a);
          boundary_flux_type &bf = d_boundary_flux[side][g][angle];
          int n1 = bf.size();
          int n2 = bf[0].size();
          for (int i = 0; i < n1; i++)
          {
            for (int j = 0; j < n2; j++)
            {
              bf[i][j] = *v;
              v++;
            }
          }
        }
      }
    }
  }
}

template <>
inline void BoundarySN<_2D>::psi(const size_t g,
                                 double *v,
                                 const int inout,
                                 const int gs,
                                 bool onlyref)
{
  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 4; side++)
  {
    if ( (not onlyref) or (onlyref and d_is_reflective[side]) )
    {
      for (int o = 0; o < 2; o++)
      {
        int octant = d_quadrature->incident_octant(side)[o];
        if (inout == OUT) octant = d_quadrature->outgoing_octant(side)[o];

        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(octant, a);
          boundary_flux_type &bf = d_boundary_flux[side][g][angle];
          int n1 = bf.size();
          for (int i = 0; i < n1; i++)
          {
            bf[i] = *v;
            v++;
          }
        }
      }
    }
  }
}

template <>
inline void BoundarySN<_1D>::psi(const size_t g,
                                 double *v,
                                 const int inout,
                                 const int gs,
                                 bool onlyref)
{

  // Incident octants for each side.
  int io[2][1] = {0, 1};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 2; side++)
  {
    if ( (not onlyref) or (onlyref and d_is_reflective[side]) )
    {
      for (int o = 0; o < 1; o++)
      {
        int octant = d_quadrature->incident_octant(side)[o];
        if (inout == OUT) octant = d_quadrature->outgoing_octant(side)[o];

        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {

          int angle = d_quadrature->index(octant, a);
          d_boundary_flux[side][g][angle] = *v;
          v++;
        }
      }
    }
  }
}

//---------------------------------------------------------------------------//
// Instantiations
//---------------------------------------------------------------------------//

template class BoundarySN<_1D>;
template class BoundarySN<_2D>;
template class BoundarySN<_3D>;

} // end namespace detran

#endif /* BOUNDARYSN_I_HH_ */
