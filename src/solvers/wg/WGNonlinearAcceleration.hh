//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  WGNonlinearAcceleration.hh
 *  @brief WGNonlinearAcceleration class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_WGNONLINEARACCELERATION_HH_
#define detran_WGNONLINEARACCELERATION_HH_

#include "boundary/BoundaryBase.hh"
#include "callow/LinearSolverCreator.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/FissionSource.hh"
#include "transport/BoundaryTally.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"
#include "solvers/solvers_export.hh"

namespace detran
{

/**
 *  @class WGNonlinearAcceleration
 *  @brief Base class for within-group acceleration schemes.
 */
template <class D>
class WGNonlinearAcceleration
{

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  typedef detran_angle::Quadrature::SP_quadrature     SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_externalsource;
  typedef detran_external_source::
          ExternalSource::vec_externalsource          vec_externalsource;
  typedef FissionSource::SP_fissionsource             SP_fissionsource;
  typedef detran_utilities::size_t                    size_t;
  typedef typename Sweeper<D>::SP_sweeper             SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource     SP_sweepsource;
  typedef State::moments_type                         moments_type;
  typedef BoundaryTally<D>::SP_tally                  SP_tally;
  typedef callow::LinearSolver::SP_solver             SP_solver;
  typedef callow::MatrixBase::SP_matrix               SP_matrix;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param external_source   User-defined external source.
   *  @param fission_source    Fission source.
   *  @param multiply          Flag for fixed source multiplying problem
   */
  WGNonlinearAcceleration(SP_state                  state,
                          SP_material               material,
                          SP_boundary               boundary,
                          SP_sweeper                sweeper,
                          SP_sweepsource            sweepsource,
                          bool                      multiply)
    : d_state(state)
    , d_material(material)
    , d_boundary(boundary)
    , d_sweeper(sweeper)
    , d_sweepsource(sweepsource)
    , d_multiply(multiply)
  {
    Require(d_state);
    Require(d_material);
    Require(d_boundary);
    Require(d_sweeper);
    Require(d_sweepsource);
    Require(d_multiply);
    d_input = d_state->get_input();
    d_mesh  = d_state->get_mesh();
  }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL WITHIN-GROUP ACCELERATORS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /**
   *  @brief Update the within-group flux moments
   *  @param  g     Group of current withing-group problem
   *  @param  phi   Flux moments to update
   */
  void update(const size_t g, moments_type &phi) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_state        d_state;
  SP_material     d_material;
  SP_boundary     d_boundary;
  SP_sweeper      d_sweeper;
  SP_sweepsource  d_sweepsource;
  bool            d_multiply;
  SP_input        d_input;
  SP_mesh         d_mesh;
  SP_tally        d_tally;
  SP_solver       d_solver;

};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//


#endif /* detran_WGNONLINEARACCELERATION_HH_ */
