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
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/FissionSource.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"

namespace detran
{

/**
 *  @class WGNonlinearAcceleration
 *  @brief Base class for within-group acceleration schemes
 */
template <class D>
class WGNonlinearAcceleration
{

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

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


  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param mat               Material definitions.
   *  @param quadrature        Angular mesh.
   *  @param boundary          Boundary fluxes.
   *  @param external_source   User-defined external source.
   *  @param fission_source    Fission source.
   *  @param multiply          Flag for fixed source multiplying problem
   */
  WGNonlinearAcceleration(SP_material               material,
                          SP_quadrature             quadrature,
                          SP_boundary               boundary,
                          const vec_externalsource &q_e,
                          SP_fissionsource          q_f,
                          bool                      multiply);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL WITHIN-GROUP SOLVERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /**
   *  @brief Update the within-group flux moments
   *  @param  g     Group of current withing-group problem
   *  @param  phi   Flux moments to update
   */
  void update(const size_t g, moments_type &phi);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_state    d_state;
  SP_material d_material;
};

} // namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//


#endif /* detran_WGNONLINEARACCELERATION_HH_ */
