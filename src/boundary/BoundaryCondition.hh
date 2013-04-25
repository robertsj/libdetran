//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryCondition.hh
 *  @author robertsj
 *  @date   Apr 9, 2012
 *  @brief  BoundaryCondition class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYCONDITION_HH_
#define detran_BOUNDARYCONDITION_HH_

#include "BoundarySN.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class BoundarySN;
template <class D> struct BoundaryTraits;

//---------------------------------------------------------------------------//
/**
 *  @class BoundaryCondition
 *  @brief Boundary condition for a surface.
 */
//---------------------------------------------------------------------------//

template <class D>
class BoundaryCondition
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<BoundaryCondition>   SP_bc;
  typedef BoundarySN<D>                             Boundary_T;
  typedef typename Boundary_T::SP_boundary          SP_boundary;
  typedef typename Boundary_T::SP_input             SP_input;
  typedef typename Boundary_T::SP_mesh              SP_mesh;
  typedef typename Boundary_T::SP_quadrature        SP_quadrature;
  typedef typename BoundaryTraits<D>::value_type    boundary_flux_type;
  typedef typename std::vector<boundary_flux_type>  vec_boundary_flux;
  typedef typename std::vector<vec_boundary_flux>   vec2_boundary_flux;
  typedef typename std::vector<vec2_boundary_flux>  vec3_boundary_flux;
  typedef detran_utilities::size_t                  size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  BoundaryCondition(SP_boundary boundary,
                    const size_t side,
                    SP_input input,
                    SP_mesh mesh,
                    SP_quadrature quadrature)
    : d_boundary(boundary)
    , d_side(side)
    , d_input(input)
    , d_mesh(mesh)
    , d_quadrature(quadrature)
  {
    Require(d_side <= D::dimension*2);
    Require(d_mesh);
    Require(d_quadrature);
  }

  /// Virtual destructor
  virtual ~BoundaryCondition(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition.
  virtual void set(const size_t g) = 0;

  /// Update a boundary following a sweep.
  virtual void update(const size_t g) = 0;

  /// Update a boundary for a given angle following a sweep.
  virtual void update(const size_t g, const size_t o, const size_t a) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Clear any static data
  virtual void clear()
  {
    /* ... */
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Boundary flux container.
  SP_boundary d_boundary;
  /// My surface.
  const size_t d_side;
  /// Input
  SP_input d_input;
  /// Cartesian mesh.
  SP_mesh d_mesh;
  /// Angular quadrature.
  SP_quadrature d_quadrature;

};

} // end namespace detran

#endif /* detran_BOUNDARYCONDITION_HH_ */
