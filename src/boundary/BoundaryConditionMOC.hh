//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryConditionMOC.hh
 * \brief  BoundaryConditionMOC 
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYCONDITIONMOC_HH_
#define BOUNDARYCONDITIONMOC_HH_

#include "BoundaryMOC.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class BoundaryMOC;

//---------------------------------------------------------------------------//
/*!
 * \class BoundaryCondition
 * \brief Boundary condition for a surface.
 */
//---------------------------------------------------------------------------//

/// \todo template on the Boundary type!!

template <class D>
class BoundaryConditionMOC
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<BoundaryConditionMOC>    SP_bc;
  typedef BoundaryMOC<D>                                Boundary_T;
  typedef typename Boundary_T::SP_boundary              SP_boundary;
  typedef typename Boundary_T::SP_input                 SP_input;
  typedef typename Boundary_T::SP_mesh                  SP_mesh;
  typedef typename Boundary_T::SP_quadrature            SP_quadrature;
  typedef detran_utilities::vec_dbl                     boundary_flux_type;
  typedef typename Boundary_T::size_t                   size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  BoundaryConditionMOC(BoundaryMOC<D>& boundary,
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
  virtual ~BoundaryConditionMOC(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MOC BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition.
  virtual void set(const size_t g) = 0;

  /// Update a boundary following a sweep.
  virtual void update(const size_t g) = 0;

  /// Update a boundary for a given angle following a sweep.
  virtual void update(const size_t g, const size_t o, const size_t a) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Boundary flux container. \todo Is there a way around using a reference?
  Boundary_T& d_boundary;

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

#endif // BOUNDARYCONDITIONMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryConditionMOC.hh
//---------------------------------------------------------------------------//
