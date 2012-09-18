//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryDiffusion.i.hh
 * \author robertsj
 * \date   Sep 11, 2012
 * \brief  BoundaryDiffusion inline member definitions
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYDIFFUSION_I_HH_
#define BOUNDARYDIFFUSION_I_HH_

#include <iostream>

namespace detran
{


//---------------------------------------------------------------------------//
// Access via cardinal indices (i.e. for sweeping)
//---------------------------------------------------------------------------//

template <class D>
inline const typename BoundaryDiffusion<D>::boundary_flux_type&
BoundaryDiffusion<D>::operator()(const size_t side,
                                 const size_t g,
                                 const size_t inout) const
{
  Require(side < D::dimension * 2);
  Require(g < d_number_groups);
  Require(inout < 2);
  return d_boundary_flux[inout][side][g];
}

template <class D>
inline typename BoundaryDiffusion<D>::boundary_flux_type&
BoundaryDiffusion<D>::operator()(const size_t side,
                                 const size_t g,
                                 const size_t inout)
{
  // Cast away return type
  return const_cast<typename BoundaryDiffusion<D>::boundary_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const BoundaryDiffusion<D>&>(*this)(side, g, inout)
  );
}

template <class D>
inline void BoundaryDiffusion<D>::clear(const size_t g)
{
  THROW("NOT IMPLEMENTED");
}

template <>
inline void BoundaryDiffusion<_3D>::clear(const size_t g)
{

  // Clear the boundary flux.
  for (int inout = 0; inout < 2; inout++)
  {
    for (int side = 0; side < 6; side++)
    {
      for (int j = 0; j < d_boundary_flux[inout][side][g].size(); j++)
      {
        for (int i = 0; i < d_boundary_flux[inout][side][g][0].size(); i++)
        {
          d_boundary_flux[inout][side][g][j][i] = 0.0;
        } // end i loop
      } // end j loop
    } // end side loop
  } // end inout loop
}

template <>
inline void BoundaryDiffusion<_2D>::clear(const size_t g)
{
  // Clear the boundary flux.
  for (int inout = 0; inout < 2; inout++)
  {
    for (int side = 0; side < 4; side++)
    {
      for (int i = 0; i < d_boundary_flux[inout][side][g].size(); i++)
      {
        d_boundary_flux[inout][side][g][i] = 0.0;
      }
    }
  }
}

template <>
inline void BoundaryDiffusion<_1D>::clear(const size_t g)
{
  // Clear the boundary flux.
  for (int inout = 0; inout < 2; inout++)
  {
    for (int side = 0; side < 2; side++)
    {
      d_boundary_flux[inout][side][g] = 0.0;
    }
  }
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class BoundaryDiffusion<_1D>;
template class BoundaryDiffusion<_2D>;
template class BoundaryDiffusion<_3D>;

} // end namespace detran


#endif /* BOUNDARYDIFFUSION_I_HH_ */
