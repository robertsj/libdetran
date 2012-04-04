//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerIteration.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef INNERITERATION_HH_
#define INNERITERATION_HH_

// Detran
#include "Boundary.hh"
#include "State.hh"
#include "Sweeper2D.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class InnerIteration
 * \brief 
 */
//===========================================================================//

template <class D>
class InnerIteration
{

public:

  typedef detran_utils::SP<InnerIteration>      SP_inner;
  // basic objects
  typedef detran_utils::InputDB::SP_input       SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename Boundary<D>::SP_boundary     SP_boundary;
  // source typedefs

  //
  typedef typename SweeperBase<D>::SP_sweeper   SP_sweeper;

  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param external_source   User-defined external source.
   *  \param fission_source    Fission source.
   */
  InnerIteration(SP_input       input,
                 SP_state       state,
                 SP_mesh        mesh,
                 SP_material    material,
                 SP_quadrature  quadrature,
                 SP_boundary    boundary);

  /*!
   *  \brief Solve the within group equation.
   */
  virtual void solve(int g) = 0;


private:

  /// User input.
  SP_input d_input;
  /// State vectors.
  SP_state d_state;
  /// Problem mesh (either Cartesian mesh or MOC tracking)
  SP_mesh d_mesh;
  /// Materials definitions.
  SP_material d_material;
  /// Angular mesh.
  SP_quarature d_quadrature;
  /// Boundary fluxes.
  SP_boundary d_boundary;


  /// Sweeper over the space-angle domain.
  SP_sweeper d_sweeper;
  /// Sweep source
  //d_sweep_source
  /// User-defined external source
  //SP_externalsource d_external_source;
  /// Fission source, if used
  //SP_fissionsource d_fission_source;
  /// Any source that remains "fixed" within the group solve
  //d_fixed_source
  /// The within group scattering source
  //SP_scattersource d_scatter_source;
  /// Maximum iterations
  int d_max_iters;
  /// Convergence tolerance
  double d_tolerance;
  /// Group we are solving.
  int d_g
  /// Moments to discrete operator.
  //d_M
  /// Numer of sweeps
  int d_number_sweeps;
};

} // end namespace detran

#endif /* INNERITERATION_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerIteration.hh
//---------------------------------------------------------------------------//
