//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Reflective.hh
 * \author robertsj
 * \date   Apr 9, 2012
 * \brief  Reflective class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef REFLECTIVE_HH_
#define REFLECTIVE_HH_

// Detran
#include "BoundaryCondition.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class Boundary;
template <class D> class BoundaryTraits;

//---------------------------------------------------------------------------//
/*!
 * \class Reflective
 * \brief Reflective boundary condition.
 */
//---------------------------------------------------------------------------//

template <class D>
class Reflective : public BoundaryCondition<D>
{

public:

  typedef SP<Reflective>                            SP_bc;
  typedef typename BoundaryCondition<D>::SP_bc      SP_base;
  typedef typename Boundary<D>::SP_boundary         SP_boundary;
  typedef InputDB::SP_input                         SP_input;
  typedef Mesh::SP_mesh                             SP_mesh;
  typedef Quadrature::SP_quadrature                 SP_quadrature;
  typedef typename BoundaryTraits<D>::value_type    boundary_flux_type;
  typedef std::vector<boundary_flux_type>           vec_boundary_flux;
  typedef std::vector<vec_boundary_flux>            vec2_boundary_flux;
  typedef std::vector<vec2_boundary_flux>           vec3_boundary_flux;


  Reflective(SP_boundary boundary,
             int side,
             SP_input input,
             SP_mesh mesh,
             SP_quadrature quadrature)
    : BoundaryCondition<D>(boundary, side, input, mesh, quadrature)
  {
    /* ... */
  }

  /// Set initial and/or fixed boundary condition.
  void set(int g){}

  /// Update a boundary following a sweep.
  void update(int g){}

private:

};

template class Reflective<_1D>;
template class Reflective<_2D>;
template class Reflective<_3D>;

} // end namespace detran

#endif /* REFLECTIVE_HH_ */
