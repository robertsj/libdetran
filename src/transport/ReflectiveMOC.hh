//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ReflectiveMOC.hh
 * \brief  ReflectiveMOC class definition.
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef REFLECTIVEMOC_HH_
#define REFLECTIVEMOC_HH_

// Detran
#include "BoundaryConditionMOC.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class BoundaryMOC;

//---------------------------------------------------------------------------//
/*!
 * \class ReflectiveMOC
 * \brief Vacuum boundary condition for MOC
 */
//---------------------------------------------------------------------------//

template <class D>
class ReflectiveMOC : public BoundaryConditionMOC<D>
{

public:

  typedef BoundaryConditionMOC<D>         Base;
  typedef typename Base::SP_bc            SP_base;
  typedef InputDB::SP_input               SP_input;
  typedef SP<MeshMOC>                     SP_mesh;
  typedef SP<QuadratureMOC>               SP_quadrature;
  typedef BoundaryMOC<D>                  Boundary_T;

  ReflectiveMOC(Boundary_T& boundary,
         int side,
         SP_input input,
         SP_mesh mesh,
         SP_quadrature quadrature)
    : Base(boundary, side, input, mesh, quadrature)
  {
    /* ... */
  }

  /// Set initial and/or fixed boundary condition.  Reflective does nothing.
  void set(int g){};

  /// Update a boundary following a sweep.
  void update(int g);

  /// Update a boundary for a given angle following a sweep.
  void update(int g, int o, int a);

private:

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

#include "ReflectiveMOC.i.hh"

#endif // REFLECTIVEMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file ReflectiveMOC.hh
//---------------------------------------------------------------------------//
