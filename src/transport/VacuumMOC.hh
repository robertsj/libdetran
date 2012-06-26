//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   VacuumMOC.hh
 * \brief  VacuumMOC class definition.
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

#ifndef VACUUMMOC_HH_
#define VACUUMMOC_HH_

// Detran
#include "BoundaryConditionMOC.hh"

namespace detran
{

// Forward declare boundary and traits.
template <class D> class BoundaryMOC;

//---------------------------------------------------------------------------//
/*!
 * \class VacuumMOC
 * \brief Vacuum boundary condition for MOC
 */
//---------------------------------------------------------------------------//

template <class D>
class VacuumMOC : public BoundaryConditionMOC<D>
{

public:

  typedef InputDB::SP_input   SP_input;
  typedef SP<MeshMOC>         SP_mesh;
  typedef SP<QuadratureMOC>   SP_quadrature;

  VacuumMOC(BoundaryMOC<D>& boundary,
         int side,
         SP_input input,
         SP_mesh mesh,
         SP_quadrature quadrature)
    : BoundaryConditionMOC<D>(boundary, side, input, mesh, quadrature)
  {
    /* ... */
  }

  /// Set initial and/or fixed boundary condition.  Vacuum does nothing.
  void set(int g){}

  /// Update a boundary following a sweep.  Vacuum does nothing.
  void update(int g){}

  ///
  void update(int g, int o, int a){}

private:


};


} // end namespace detran

#endif // VACUUMMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file VacuumMOC.hh
//---------------------------------------------------------------------------//
