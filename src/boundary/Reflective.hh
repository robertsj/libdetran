//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Reflective.hh
 *  @author robertsj
 *  @date   Apr 9, 2012
 *  @brief  Reflective class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_REFLECTIVE_HH_
#define detran_REFLECTIVE_HH_

#include "BoundaryCondition.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Reflective
 *  @brief Reflective boundary condition so SN problems.
 */
//---------------------------------------------------------------------------//

template <class D>
class Reflective : public BoundaryCondition<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryCondition<D>                Base;
  typedef typename Base::SP_bc                SP_bc;
  typedef typename Base::Boundary_T           Boundary_T;
  typedef typename Base::SP_boundary          SP_boundary;
  typedef typename Base::SP_input             SP_input;
  typedef typename Base::SP_mesh              SP_mesh;
  typedef typename Base::SP_quadrature        SP_quadrature;
  typedef typename Base::vec_boundary_flux    vec_boundary_flux;
  typedef typename Base::vec2_boundary_flux   vec2_boundary_flux;
  typedef typename Base::vec3_boundary_flux   vec3_boundary_flux;
  typedef typename Base::size_t               size_t;
  typedef detran_utilities::vec_int           vec_int;
  typedef detran_utilities::vec2_int          vec2_int;
  typedef detran_geometry::Mesh               Mesh;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Reflective(SP_boundary boundary,
             const size_t side,
             SP_input input,
             SP_mesh mesh,
             SP_quadrature quadrature)
    : Base(boundary, side, input, mesh, quadrature)
    , d_octants(quadrature->number_octants()/2, vec_int(2, 0))
  {
    setup_octant();
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition.  Reflective does nothing.
  void set(const size_t g){};

  /// Update a boundary following a sweep.
  void update(const size_t g);

  /// Update a boundary for a given angle following a sweep.
  void update(const size_t g, const size_t o, const size_t a);

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Index of octant reflection pairs.
  vec2_int d_octants;

  // Set up the octant indices.
  void setup_octant();

  // Make inherited data visible
  using Base::d_boundary;
  using Base::d_side;
  using Base::d_input;
  using Base::d_mesh;
  using Base::d_quadrature;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Reflective.i.hh"

#endif /* detran_REFLECTIVE_HH_ */
