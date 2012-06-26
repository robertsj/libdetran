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

// Detran
#include "BoundaryMOC.hh"
#include "QuadratureMOC.hh"
#include "MeshMOC.hh"
#include "Traits.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "InputDB.hh"
#include "SP.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class BoundaryMOC;
//template <class D> class BoundaryTraits;

//---------------------------------------------------------------------------//
/*!
 * \class BoundaryCondition
 * \brief Boundary condition for a surface.
 */
//---------------------------------------------------------------------------//

/// \todo template on the Boundary type!!

template <class D>
class BoundaryConditionMOC : public Object
{

public:

  typedef SP<BoundaryConditionMOC>                  SP_bc;
  typedef typename BoundaryMOC<D>::SP_boundary      SP_boundary;
  typedef InputDB::SP_input                         SP_input;
  typedef MeshMOC::SP_mesh                          SP_mesh;
  typedef QuadratureMOC::SP_quadrature              SP_quadrature;
  typedef vec_dbl                                   boundary_flux_type;

  BoundaryConditionMOC(BoundaryMOC<D>& boundary,
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

  /// Virtual destructor
  virtual ~BoundaryConditionMOC(){}

  /// Set initial and/or fixed boundary condition.
  virtual void set(int g) = 0;

  /// Update a boundary following a sweep.
  virtual void update(int g) = 0;

  /// Update a boundary for a given angle following a sweep.
  virtual void update(int g, int o, int a) = 0;

  /// DBC function
  virtual bool is_valid() const
  {
    return true;
  }

protected:

  /// Boundary flux container. \todo Is there a way around using a reference?
  BoundaryMOC<D>& d_boundary;

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

#endif // BOUNDARYCONDITIONMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryConditionMOC.hh
//---------------------------------------------------------------------------//
