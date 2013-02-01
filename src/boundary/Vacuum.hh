//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Vacuum.hh
 *  @author robertsj
 *  @date   Apr 9, 2012
 *  @brief  Vacuum class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_VACUUM_HH_
#define detran_VACUUM_HH_

#include "BoundaryCondition.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Vacuum
 *  @brief Vacuum boundary condition for SN calculations
 */
//---------------------------------------------------------------------------//

template <class D>
class Vacuum : public BoundaryCondition<D>
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
  typedef typename Base::size_t               size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Vacuum(SP_boundary boundary,
         const size_t side,
         SP_input input,
         SP_mesh mesh,
         SP_quadrature quadrature)
    : Base(boundary, side, input, mesh, quadrature)
  {
    /* ... */
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition. Vacuum does nothing.
  void set(const size_t g){}

  /// Update a boundary following a sweep. Vacuum does nothing.
  void update(const size_t g){}

  /// Update a boundary for a given angle following a sweep. Vacuum does nothing.
  void update(const size_t g, const size_t o, const size_t a){}

private:

};

} // end namespace detran

#endif /* detran_VACUUM_HH_ */
