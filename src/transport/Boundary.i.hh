//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Boundary.i.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Boundary inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef BOUNDARY_I_HH_
#define BOUNDARY_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
// Access via cardinal indices (i.e. for sweeping)
//---------------------------------------------------------------------------//

template <class D>
inline const typename Boundary<D>::boundary_flux_type&
Boundary<D>::operator()(int side, int o, int a, int g) const
{
  Require(side >= 0);
  Require(side < D::dimension*2);
  Require(d_quadrature->valid_index(o, a));
  Require(g >= 0);
  Require(g < d_number_groups);
  int angle = d_quadrature->index(o, a);
  return d_boundary_flux[side][g][angle];
}

template <class D>
inline typename Boundary<D>::boundary_flux_type&
Boundary<D>::operator()(int side, int o, int a, int g)
{
  // Cast away return type
  return const_cast<typename Boundary<D>::boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const Boundary<D>&>(*this)(side, o, a, g)
  );
}

//---------------------------------------------------------------------------//
// Access to incident and outgoing flux with local ordering
//---------------------------------------------------------------------------//

template <class D>
inline const typename Boundary<D>::boundary_flux_type&
Boundary<D>::incident(int side, int angle, int g) const
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
inline typename Boundary<D>::boundary_flux_type&
Boundary<D>::incident(int side, int angle, int g)
{
  // Cast away return type
  return const_cast<typename Boundary<D>::boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const Boundary<D>&>(*this).incident(side, angle, g)
  );
}

template <class D>
inline const typename Boundary<D>::boundary_flux_type&
Boundary<D>::outgoing(int side, int angle, int g) const
{
  int index = ordered_angle(side, angle, OUT);
}

template <class D>
inline typename Boundary<D>::boundary_flux_type&
Boundary<D>::outgoing(int side, int angle, int g)
{
  // Cast away return type
  return const_cast<boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const Boundary<D>&>(*this).outgoing(side, angle, g)
  );
}

template <class D>
inline void Boundary<D>::set(int g)
{
  // Set boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->set(g);
  }
}

template <class D>
inline void Boundary<D>::update(int g)
{
  // Set boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->update(g);
  }
}



//---------------------------------------------------------------------------//
// Local ordering index.
//---------------------------------------------------------------------------//

template <class D>
inline int Boundary<D>::ordered_angle(int side, int angle, int inout) const
{
  return 0;
}

template <>
inline int Boundary<_2D>::ordered_angle(int side, int angle, int inout) const
{
  return 0;
}

template <>
inline int Boundary<_1D>::ordered_angle(int side, int angle, int inout) const
{
  //
  if (inout == IN)
  {

  }

}

//---------------------------------------------------------------------------//
// Instantiations
//---------------------------------------------------------------------------//

template class BoundaryTraits<_1D>;
template class BoundaryTraits<_2D>;
template class BoundaryTraits<_3D>;

template class Boundary<_1D>;
template class Boundary<_2D>;
template class Boundary<_3D>;

} // end namespace detran

#endif /* BOUNDARY_I_HH_ */
