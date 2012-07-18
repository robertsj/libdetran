//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryCondition.hh
 * \author robertsj
 * \date   Apr 9, 2012
 * \brief  BoundaryCondition class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYCONDITION_HH_
#define BOUNDARYCONDITION_HH_

// Detran
#include "Boundary.hh"
#include "Quadrature.hh"
#include "Mesh.hh"
#include "Traits.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "InputDB.hh"
#include "SP.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class Boundary;
template <class D> class BoundaryTraits;

//---------------------------------------------------------------------------//
/*!
 * \class BoundaryCondition
 * \brief Boundary condition for a surface.
 */
//---------------------------------------------------------------------------//

template <class D>
class BoundaryCondition : public Object
{

public:

  typedef SP<BoundaryCondition>                     SP_bc;
  typedef Boundary<D>                               Boundary_T;
  typedef typename Boundary_T::SP_boundary          SP_boundary;
  typedef InputDB::SP_input                         SP_input;
  typedef Mesh::SP_mesh                             SP_mesh;
  typedef Quadrature::SP_quadrature                 SP_quadrature;
  typedef typename BoundaryTraits<D>::value_type    boundary_flux_type;
  typedef std::vector<boundary_flux_type>           vec_boundary_flux;
  typedef std::vector<vec_boundary_flux>            vec2_boundary_flux;
  typedef std::vector<vec2_boundary_flux>           vec3_boundary_flux;


  BoundaryCondition(Boundary_T& boundary,
                    int side,
                    SP_input input,
                    SP_mesh mesh,
                    SP_quadrature quadrature)
    : d_boundary(boundary)
    , d_side(side)
    , d_input(input)
    , d_mesh(mesh)
    , d_quadrature(quadrature)
  {
    //Require(d_boundary);
    Require(d_side >= 0);
    Require(d_side <= D::dimension*2);
    Require(d_mesh);
    Require(d_quadrature);
  }

  /// Virtual destructor so Eclipse stops complaining.
  virtual ~BoundaryCondition()
  {

  }

  /// Set initial and/or fixed boundary condition.
  virtual void set(int g) = 0;

  /// Update a boundary following a sweep.
  virtual void update(int g) = 0;

  /// Update a boundary for a given angle following a sweep.
  virtual void update(int g, int o, int a) = 0;

  virtual bool is_valid() const
  {return true;}

protected:

  /// Boundary flux container. \todo Is there a way around using a reference?
  Boundary_T& d_boundary;

  /// My surface.
  const int d_side;

  /// Input
  SP_input d_input;

  /// Cartesian mesh.
  SP_mesh d_mesh;

  /// Angular quadrature.
  SP_quadrature d_quadrature;
};

} // end namespace detran

#endif /* BOUNDARYCONDITION_HH_ */
