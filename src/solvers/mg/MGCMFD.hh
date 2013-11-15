//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCMFD.hh
 *  @brief MGCMFD class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGCMFD_HH_
#define detran_MGCMFD_HH_

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class MGCMFD
 *  @brief Creates operators for multigroup CMFD fixed and eigenvalue problems
 *
 *  MGCMFD constructs the coarse mesh and current tally required to build
 *  CMFD diffusion operators.  These operators can be used
 */
//----------------------------------------------------------------------------//
/**
 *  @example solvers/test/test_MGCMFD.cc
 *
 *  Test of MGCMFD class.
 */
template <class D>
class MGCMFD
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef State::SP_state                           SP_state;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary     SP_boundary;
  typedef typename Sweeper<D>::SP_sweeper           SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource   SP_sweepsource;
  typedef detran_utilities::SP<callow::Matrix>      SP_operator;
  typedef callow::Vector::SP_vector                 SP_vector;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::size_t                  size_t;
  typedef detran_utilities::vec_size_t              groups_t;
  typedef groups_t::iterator                        groups_iter;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param q_e               Vector of user-defined external sources
   *  @param q_f               Fission source.
   *  @param multiply          Flag for a multiplying fixed source problem
   */
  MGCMFD(SP_state                   state,
         SP_material                material,
         SP_boundary                boundary,
         SP_sweeper                 sweeper,
         SP_sweepsource             sweepsource,
         bool                       multiply = false);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Return dtilde
  vec2

  /// Return the CMFD loss operator
  SP_operator get_loss_operator(double keff = 1.0);

  /// Return the CMFD gain operator
  SP_operator get_gain_operator();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// State
  SP_state d_state;
  /// Material
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// Boundary
  SP_boundary d_boundary;
  /// Sweeper
  SP_sweeper d_sweeper;
  /// Sweep source
  SP_sweepsource d_sweepsource;
  ///
  bool d_multiply;
  size_t d_dimension;
  size_t d_number_groups;
  size_t d_group_size;
  /// Loss operator
  SP_operator d_loss;
  /// Gain operator
  SP_operator d_gain;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Build the right hand side.
  void build_rhs(State::vec_moments_type &B);

  /**
   *  @brief Sweep through the upscatter block once
   *
   *  For computing the right hand side (i.e. the uncollided flux), we
   *  need to sweep through all the groups once.  Within a group, we need
   *  to iterate on reflecting conditions.  This process is also need at
   *  the end of a solve to pick up the boundary fluxes.
   *
   *  @param phi    multigroup fluxes to update
   */
  void group_sweep(State::vec_moments_type &phi);

  /// Solve the within-group reflection problem
  void solve_wg_reflection(size_t g, State::moments_type &phi_g);

};

} // namespace detran

#endif /* detran_MGCMFD_HH_ */

//----------------------------------------------------------------------------//
//              end of MGCMFD.hh
//----------------------------------------------------------------------------//
