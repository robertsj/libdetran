//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vacuum.hh
 * \author robertsj
 * \date   Apr 9, 2012
 * \brief  Vacuum class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef VACUUM_HH_
#define VACUUM_HH_

// Detran
#include "BoundaryCondition.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class Boundary;
template <class D> class BoundaryTraits;

//---------------------------------------------------------------------------//
/*!
 * \class Vacuum
 * \brief Vacuum boundary condition.
 */
//---------------------------------------------------------------------------//

template <class D>
class Vacuum : public BoundaryCondition<D>
{

public:

  typedef SP<Vacuum>                                SP_bc;
  typedef typename BoundaryCondition<D>::SP_bc      SP_base;
  typedef typename Boundary<D>::SP_boundary         SP_boundary;
  typedef InputDB::SP_input                         SP_input;
  typedef Mesh::SP_mesh                             SP_mesh;
  typedef Quadrature::SP_quadrature                 SP_quadrature;
  typedef typename BoundaryTraits<D>::value_type    boundary_flux_type;
  typedef std::vector<boundary_flux_type>           vec_boundary_flux;
  typedef std::vector<vec_boundary_flux>            vec2_boundary_flux;
  typedef std::vector<vec2_boundary_flux>           vec3_boundary_flux;


  Vacuum(Boundary<D>& boundary,
         int side,
         SP_input input,
         SP_mesh mesh,
         SP_quadrature quadrature)
    : BoundaryCondition<D>(boundary, side, input, mesh, quadrature)
  {
    /* ... */
  }

  /// Set initial and/or fixed boundary condition.  Vacuum does nothing.
  void set(int g){}

  /// Update a boundary following a sweep.  Vacuum does nothing.
  void update(int g){}

  ///
  void update(int g, int o, int a){}

private:

};

//template class Vacuum<_1D>;
//template class Vacuum<_2D>;
//template class Vacuum<_3D>;

} // end namespace detran

#endif /* VACUUM_HH_ */
