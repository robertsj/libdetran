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

#include "BoundaryConditionMOC.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class ReflectiveMOC
 * \brief Reflective boundary condition for MOC
 */
//---------------------------------------------------------------------------//

template <class D>
class ReflectiveMOC : public BoundaryConditionMOC<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryConditionMOC<D>             Base;
  typedef typename Base::SP_bc                SP_bc;
  typedef typename Base::Boundary_T           Boundary_T;
  typedef typename Base::SP_input             SP_input;
  typedef typename Base::SP_mesh              SP_mesh;
  typedef typename Base::SP_quadrature        SP_quadrature;
  typedef detran_utilities::vec_int           vec_int;
  typedef detran_utilities::vec2_int          vec2_int;
  typedef typename Base::size_t               size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  ReflectiveMOC(Boundary_T& boundary,
                const size_t side,
                SP_input input,
                SP_mesh mesh,
                SP_quadrature quadrature)
    : Base(boundary, side, input, mesh, quadrature)
  {
    /* ... */
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MOC BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition.  Reflective does nothing.
  void set(const size_t g){};

  /// Update a boundary following a sweep.
  void update(const size_t g);

  /// Update a boundary for a given angle following a sweep.
  void update(const size_t g, const size_t o, const size_t a);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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
