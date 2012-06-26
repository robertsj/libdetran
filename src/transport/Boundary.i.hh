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

// System
#include <iostream>

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

template <class D>
inline void Boundary<D>::update(int g, int o, int a)
{
  // Set boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    d_bc[side]->update(g, o, a);
  }
}

template <class D>
inline void Boundary<D>::clear(int g)
{
  // Clear the boundary flux.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    for (int angle = 0; angle < d_boundary_flux[side][g].size(); angle++)
    {
      for (int j = 0; j < d_mesh->number_cells_y(); j++)
      {
        for (int i = 0; i < d_mesh->number_cells_x(); i++)
        {
          d_boundary_flux[side][g][angle][i][j] = 0.0;
        }
      }
    }
  }
}

template <>
inline void Boundary<_2D>::clear(int g)
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
inline void Boundary<_1D>::clear(int g)
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

// Incident condition setter and getter for Krylov solvers.

template<class D>
inline void Boundary<D>::set_incident(int g, double *v)
{
  /* ... */
}

template<>
inline void Boundary<_3D>::set_incident(int g, double *v)
{
  // Incident octants for each side.
  int io[6][4] = {0,3,4,7,  2,1,5,6,  0,4,1,5,  7,3,6,2,  1,0,2,3,  4,5,7,6};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 6; side++)
  {
    if (d_is_reflective[side])
    {
      for (int o = 0; o < 4; o++)
      {
        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(io[side][o], a);
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

template<>
inline void Boundary<_2D>::set_incident(int g, double *v)
{
  // Incident octants for each side.
  int io[4][2] = {0,3,  2,1,  1,0,  3,2};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 4; side++)
  {
    if (d_is_reflective[side])
    {
      for (int o = 0; o < 2; o++)
      {
        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(io[side][o], a);
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

template<>
inline void Boundary<_1D>::set_incident(int g, double *v)
{
  using std::cout;
  using std::endl;

  // Incident octants for each side.
  int io[2][1] = {0, 1};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 2; side++)
  {
    if (d_is_reflective[side])
    {
      for (int o = 0; o < 1; o++)
      {
        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {

          int angle = d_quadrature->index(io[side][o], a);
          d_boundary_flux[side][g][angle] = *v;
          v++;
        }
      }
    }
  }
}

template<class D>
inline void Boundary<D>::get_incident(int g, double *v)
{
  /* ... */
}

template<>
inline void Boundary<_3D>::get_incident(int g, double *v)
{
  // Incident octants for each side.
  int io[6][4] = {0,3,4,7,  2,1,5,6,  0,4,1,5,  7,3,6,2,  1,0,2,3,  4,5,7,6};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 6; side++)
  {
    if (d_is_reflective[side])
    {
      for (int o = 0; o < 4; o++)
      {
        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(io[side][o], a);
          boundary_flux_type &bf = d_boundary_flux[side][g][angle];
          int n1 = bf.size();
          int n2 = bf[0].size();
          for (int i = 0; i < n1; i++)
          {
            for (int j = 0; j < n2; j++)
            {
              *v = bf[i][j];
              v++;
            }
          }
        }
      }
    }
  }
}

template<>
inline void Boundary<_2D>::get_incident(int g, double *v)
{
  // Incident octants for each side.
  int io[4][2] = {0,3,  2,1,  1,0,  3,2};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 4; side++)
  {
    if (d_is_reflective[side])
    {
      for (int o = 0; o < 2; o++)
      {
        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(io[side][o], a);
          boundary_flux_type &bf = d_boundary_flux[side][g][angle];
          int n1 = bf.size();
          for (int i = 0; i < n1; i++)
          {
            *v = bf[i];
            v++;
          }
        }
      }
    }
  }
}

template<>
inline void Boundary<_1D>::get_incident(int g, double *v)
{
  // Incident octants for each side.
  int io[2][1] = {0, 1};

  // Loop over all reflective sides and set the flux.
  for (int side = 0; side < 2; side++)
  {
    if (d_is_reflective[side])
    {
      for (int o = 0; o < 1; o++)
      {
        for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
        {
          int angle = d_quadrature->index(io[side][o], a);
          *v = d_boundary_flux[side][g][angle];
          v++;
        }
      }
    }
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
