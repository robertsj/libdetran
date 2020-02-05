//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SweepSource.hh
 *  @brief SweepSource class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SWEEPSOURCE_HH_
#define detran_SWEEPSOURCE_HH_

#include "DimensionTraits.hh"
#include "FissionSource.hh"
#include "ScatterSource.hh"
#include "State.hh"
#include "angle/MomentToDiscrete.hh"
#include "angle/Quadrature.hh"
#include "external_source/ExternalSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran
{

//----------------------------------------------------------------------------//
/**
 *  @class SweepSource
 *  @brief Construct the source for sweeping
 *
 *  This class defines a general right-hand-side construct for sweeps (which
 *  are the basis for all Detran transport solvers).  That is, the SweepSource
 *  is a source for a sweep along a particular angle.
 *
 *  The source for a given sweep is comprised in general of these four
 *  components:
 *    - within-group scattering source
 *    - in-scattering source (from both up- and down-scattering)
 *    - fission source
 *    - external source
 *  Algorithmically, the fission source behaves just like scattering for
 *  fixed source problems unless explicitly treated as an external source.  For
 *  eigenvalue problems, it is always treated as an external source.
 *
 *  Recall that a within group equation is represented as
 *  @f[
 *      \mathbf{T}[\psi]_g =
 *          [\mathbf{M}][\mathbf{S}]_{gg}[\phi]_g + \bar{Q}_g \, ,
 *  @f]
 *  where \f$ \bar{Q}_g \f$ represents everything \em but the within-group
 *  scattering (and fission, if applicable).
 *  The right hand side is the sweep source (but \em not the
 *  right hand side for the linear system of interest, which casts the
 *  problem in terms of flux moments!).
 *
 *  As a sanity check, consider a purely
 *  isotropic flux in one group with isotropic scattering. Then we must
 *  have the equation
 *  @f[
 *      \mathbf{T}\psi = \Sigma_{0} \phi_{00} / 4\pi \, ,
 *  @f]
 *  which indicates that application of \f$ \mathbf{M} \f$ implies
 *  a normalization.
 *
 *  @sa ScatterSource, FissionSource, ExternalSource
 */
//----------------------------------------------------------------------------//

template<class D>
class SweepSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<SweepSource>         SP_sweepsource;
  typedef State::SP_state                           SP_state;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef detran_angle::MomentToDiscrete::SP_MtoD   SP_MtoD;
  typedef detran_external_source::
          ExternalSource::SP_externalsource         SP_externalsource;
  typedef ScatterSource::SP_scattersource           SP_scattersource;
  typedef FissionSource::SP_fissionsource           SP_fissionsource;
  typedef detran_utilities::vec_dbl                 sweep_source_type;
  typedef detran_utilities::size_t                  size_t;
  typedef State::moments_type                       moments_type;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state              State vector
   *  @param mesh               Problem geometry
   *  @param quadrature         Angular mesh
   *  @param material           Cross section definitions
   *  @param MtoD               Moment-to-discrete operator
   *  @param implicit_fission   True if fission treated like scattering
   */
  SweepSource(SP_state      state,
              SP_mesh       mesh,
              SP_quadrature quadrature,
              SP_material   material,
              SP_MtoD       MtoD,
              bool          implicit_fission = false)
    :  d_state(state)
    ,  d_mesh(mesh)
    ,  d_quadrature(quadrature)
    ,  d_MtoD(MtoD)
    ,  d_source(mesh->number_cells(), 0.0)
    ,  d_fixed_group_source(mesh->number_cells(), 0.0)
    ,  d_scatter_group_source(mesh->number_cells(), 0.0)
    ,  d_implicit_fission(implicit_fission)
    ,  d_scattersource(new ScatterSource(mesh, material, state))
    ,  d_discrete_external_source_flag(true)
  {
    Require(d_state);
    Require(d_mesh);
    Require(d_quadrature);
    Require(d_MtoD);
  }

  /// SP Constructor
  static SP_sweepsource Create(SP_state      state,
                               SP_mesh       mesh,
                               SP_quadrature quadrature,
                               SP_material   material,
                               SP_MtoD       MtoD,
                               bool          multiply = false)
  {
    SP_sweepsource p(new SweepSource(state, mesh, quadrature, material,
                                     MtoD, multiply));
    return p;
  }

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Set an external moment source.
   *  @param source  Smart pointer to external source
   */
  void set_moment_source(SP_externalsource source)
  {
    Require(source);
    d_moment_external_sources.push_back(source);
  }

  /**
   *  @brief Set an external discrete source.
   *
   *  There really is no difference between a moment and discrete
   *  source in representation (the functions are the same), but
   *  the sweep source adds the discrete portion of the source
   *  before giving the client the source.  This way, the moment
   *  sources can be done once and the discrete component can
   *  be built as needed.
   *
   *  @param source  Smart pointer to external source
   */
  void set_discrete_source(SP_externalsource source)
  {
    Require(source);
    d_discrete_external_sources.push_back(source);
  }

  /**
   *  @brief Set a fission source.
   *  @param source  Smart pointer to fission source
   */
  void set_fission_source(SP_fissionsource source)
  {
    Require(source);
    d_fissionsource = source;
  }

  SP_scattersource get_scatter_source()
  {
    return d_scattersource;
  }

  SP_fissionsource get_fission_source()
  {
    return d_fissionsource;
  }


  // The following source construction routines aim to
  // satisfy the needs of all solvers.

  /**
   *  @brief Add all fixed moments source.
   *
   *  Creates a fixed group source using external moments sources
   *  and, if applicable, fission sources.  Note, fission is
   *  only fixed within an eigenproblem.
   */
  void build_fixed(const size_t g);

  /**
   *  @brief Build fixed with in-scatter.
   *
   *  Creates a fixed group source as \ref build_fixed, but adds
   *  the in-scatter source.  This is the typical construction
   *  to be called within source iteration.
   *
   */
  void build_fixed_with_scatter(const size_t g);

  /**
   *  @brief Build fixed with downs-scatter.
   *
   *  Creates a fixed group source as \ref build_fixed, but adds
   *  the in-scatter source.  This is the typical construction
   *  to be called within source iteration.
   *
   *  @param g          Group of current solve
   *  @param g_cutoff   Last group to include in downscatter source
   */
  void build_fixed_with_downscatter(const size_t g, const size_t g_cutoff);

  /**
   *  @brief Build within-group scattering source.
   *
   *  Called before each sweep.
   */
  void build_within_group_scatter(const size_t g, const moments_type &phi);

  /**
   *  @brief Build total scattering source.
   *
   *  This builds the complete group scatter source using the
   *  given multigroup flux vector.  This routine would find
   *  use in a multigroup Krylov solver.
   */
  void build_total_scatter(const size_t g,
                           const size_t g_cutoff,
                           const State::vec_moments_type &phi);

  /// Reset all the internal source vectors to zero.
  void reset();

  /// Fill a source vector
  void source(const size_t g,
              const size_t o,
              const size_t a,
              sweep_source_type& s);

  /// Return the fixed source for the current group
  const moments_type& fixed_group_source() const
  {
    return d_fixed_group_source;
  }
  moments_type& fixed_group_source()
  {
    return d_fixed_group_source;
  }

  /// Return the scatter source for the current group
  const moments_type& scatter_group_source() const
  {
    return d_scatter_group_source;
  }

  void set_discrete_external_source_flag(bool flag)
  {
    d_discrete_external_source_flag = flag;
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Problem state vector
  SP_state d_state;
  /// Geometry
  SP_mesh d_mesh;
  /// Quadrature
  SP_quadrature d_quadrature;
  /// Moment-to-discrete operator
  SP_MtoD d_MtoD;
  /// Sweep source for a given angle and group over all cells.
  sweep_source_type d_source;
  /// Fixed moments source applicable to all angles in this group.
  moments_type d_fixed_group_source;
  /// Within group scattering applicable to all angles in this group.
  moments_type d_scatter_group_source;
  /// A container of external moments sources
  std::vector<SP_externalsource> d_moment_external_sources;
  /// A container of external discrete sources
  std::vector<SP_externalsource> d_discrete_external_sources;
  /// Fission source
  SP_fissionsource d_fissionsource;
  /// Implicit fission flag
  bool d_implicit_fission;
  /// Scattering source
  SP_scattersource d_scattersource;
  /// Discrete external source flag.  Default true.
  bool d_discrete_external_source_flag;

};

} // end namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "SweepSource.i.hh"

#endif /* detran_SWEEPSOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of SweepSource.hh
//----------------------------------------------------------------------------//
