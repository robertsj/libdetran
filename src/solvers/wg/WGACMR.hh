//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  WGACMR.hh
 *  @brief WGACMR clas definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_WGACMR_HH_
#define detran_WGACMR_HH_

#include "WGNonlinearAcceleration.hh"

namespace detran
{

/**
 *  @class WGACMR
 *  @brief Within-group angular coarse mesh rebalance
 *
 *  This nonlinear acceleration scheme generalize traditional coarse
 *  mesh rebalance (CMR) by adding additional angular moments rebalance
 *  factors.  In the zeroth-order approximation, ACMR reduces to CMR.
 */
template <class D>
class WGACMR: public WGNonlinearAcceleration<D>
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef WGNonlinearAcceleration<D>          Base;
  typedef typename Base::SP_input             SP_input;
  typedef typename Base::SP_state             SP_state;
  typedef typename Base::SP_mesh              SP_mesh;
  typedef typename Base::SP_material          SP_material;
  typedef typename Base::SP_quadrature        SP_quadrature;
  typedef typename Base::SP_boundary          SP_boundary;
  typedef typename Base::SP_externalsource    SP_externalsource;
  typedef typename Base::vec_externalsource   vec_externalsource;
  typedef typename Base::SP_fissionsource     SP_fissionsource;
  typedef typename Base::size_t               size_t;
  typedef typename Base::SP_sweeper           SP_sweeper;
  typedef typename Base::SP_sweepsource       SP_sweepsource;
  typedef typename Base::moments_type         moments_type;

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
  WGACMR(SP_state                  state,
         SP_material               material,
         SP_boundary               boundary,
         const vec_externalsource &q_e,
         SP_fissionsource          q_f,
         bool                      multiply);

  /// Virtual destructor
  virtual ~WGACMR();

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL WITHIN-GROUP ACCELERATORS MUST IMPLEMENT
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

  /// Maximum angular order of the coarse mesh rebalance factors
  size_t d_order;
  /// Coarse mesh
  SP_mesh d_coarsemesh;
  /// Boundary tally
  SP_tally d_tally;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Build the ACMR operator
  void build_operator();

};

} // end namespace detran

#endif /* detran_WGACMR_HH_ */

//----------------------------------------------------------------------------//
//              end of file WGACMR.hh
//----------------------------------------------------------------------------//
