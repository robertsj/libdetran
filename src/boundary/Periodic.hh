//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Periodic.hh
 *  @brief Periodic class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//


#ifndef detran_PERIODIC_HH_
#define detran_PERIODIC_HH_

#include "BoundaryCondition.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Periodic
 *  @brief Periodic boundary condition for SN problems.
 */
//---------------------------------------------------------------------------//

template <class D>
class Periodic : public Reflective<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Reflective<D>                       Base;
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

  Periodic(SP_boundary boundary,
           const size_t side,
           SP_input input,
           SP_mesh mesh,
           SP_quadrature quadrature)
    : Base(boundary, side, input, mesh, quadrature)
  {
    if (d_side == detran_geometry::Mesh::WEST)
      d_periodic_side = detran_geometry::Mesh::EAST;
    else if (d_side == detran_geometry::Mesh::EAST)
      d_periodic_side = detran_geometry::Mesh::WEST;
    else if (d_side == detran_geometry::Mesh::SOUTH)
      d_periodic_side = detran_geometry::Mesh::NORTH;
    else if (d_side == detran_geometry::Mesh::NORTH)
      d_periodic_side = detran_geometry::Mesh::SOUTH;
    else if (d_side == detran_geometry::Mesh::BOTTOM)
      d_periodic_side = detran_geometry::Mesh::TOP;
    else
      d_periodic_side = detran_geometry::Mesh::BOTTOM;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition.  Periodic does nothing.
  void set(const size_t g){};

  /// Update a boundary following a sweep.
  void update(const size_t g);

  /// Update a boundary for a given angle following a sweep.
  void update(const size_t g, const size_t o, const size_t a);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Periodic side
  size_t d_periodic_side;

  // Make inherited data visible
  using Base::d_boundary;
  using Base::d_side;
  using Base::d_input;
  using Base::d_mesh;
  using Base::d_quadrature;
  using Base::d_octants;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Periodic.i.hh"

#endif /* detran_PERIODIC_HH_ */
