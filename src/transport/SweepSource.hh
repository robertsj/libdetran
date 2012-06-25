//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SweepSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  SweepSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPSOURCE_HH_
#define SWEEPSOURCE_HH_

// Detran
#include "ExternalSource.hh"
#include "FissionSource.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "MomentToDiscrete.hh"
#include "Quadrature.hh"
#include "ScatterSource.hh"
#include "State.hh"
#include "Traits.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

// System
#include <vector>

namespace detran
{

//===========================================================================//
/*!
 * \class SweepSource
 * \brief Construct the source for sweeping
* This class defines a general right-hand-side construct for sweeps (which
 * are the basis for all solvers).  That is, the SweepSource is a source
 * for a sweep along a particular angle.
 *
 * The source for a given sweep is comprised in general of these four
 * components:
 *  - within-group scattering source
 *  - in-scattering source (from both up- and down-scattering)
 *  - fission source
 *  - external source
 * Algorithmically, the fission source behaves as an external source within
 * an inner iteration.  Moreover, and of importance for response function
 * generation, there \em can be a fission source (i.e. multiplication)
 * together with an fixed source.
 *
 * Recall that a within group equation is represented as
 * \f[
 *    \mathbf{T}[\psi]_g = [\mathbf{M}][\mathbf{S}]_{gg}[\phi]_g + \bar{Q}_g \, ,
 * \f]
 * where \f$ \bar{Q}_g \f$ represents everything \em but the within-group
 * scattering.  The right hand side is the sweep source (but \em not the
 * right hand side for the linear system of interest, which casts the
 * problem in terms of flux moments!).
 *
 * As a sanity check, consider a purely
 * isotropic flux in one group with isotropic scattering. Then we must
 * have the equation
 * \f[
 *    \mathbf{T}\psi = \Sigma_{0} \phi_{00} / 4\pi \, ,
 * \f]
 * which indicates that application of \f$ \mathbf{M} \f$ implies
 * a normalization.
 *
 * For each inner iteration, the scattering components are stored in
 * moments form.  The M operator is then applied using one row at a
 * time for the associated angle being swept.
 *
 * \note Keep in mind where we're headed: solving the linear system
 *       \f$ (\mathbf{I}-\mathbf{D}\mathbf{T}^{-1}\mathbf{M}\mathbf{S})\phi
 *        = \mathbf{D}\mathbf{T}^{-1}q \f$. Hence, the sweep source is q,
 *        and we'll invert the transport operator \em T via sweeps.
 *
 * \note Denovo uses a different class for each source component, and then
 *       a "database" to house all the relevant components.  Here, we'll
 *       keep everything in a single class, with different methods that
 *       will set various components.  Of course, fission and external
 *       sources are still defined externally.
 *
 * \sa ScatterSource, FissionSource, ExternalSource
 */
//===========================================================================//

template<class D>
class SweepSource
{

public:

  typedef SP<SweepSource>                   SP_sweepsource;
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef typename
      MomentToDiscrete<D>::SP_MtoD          SP_MtoD;
  //
  typedef ExternalSource::SP_source         SP_externalsource;
  typedef ScatterSource::SP_source          SP_scattersource;
  typedef FissionSource::SP_source          SP_fissionsource;
  //
  typedef vec_dbl                           sweep_source_type;
  typedef State::moments_type               moments_type;

  /*!
   * \brief Constructor.
   *
   * This sizes the source variables.
   *
   * \param state           The state.
   * \param mesh            The mesh.
   * \param angularmesh     The angular mesh.
   * \param materials       The material library.
   * \param momentsindex    The moments index.
   * \param m_operator      The moments to discrete operator.
   *
   */
  SweepSource(SP_state      state,
              SP_mesh       mesh,
              SP_quadrature quadrature,
              SP_material   material,
              SP_MtoD       MtoD)
    :  d_state(state)
    ,  d_mesh(mesh)
    ,  d_quadrature(quadrature)
    ,  d_MtoD(MtoD)
    ,  d_source(mesh->number_cells(), 0.0)
    ,  d_fixed_group_source(mesh->number_cells(), 0.0)
    ,  d_scatter_group_source(mesh->number_cells(), 0.0)
    ,  d_scattersource(new ScatterSource(mesh, material, state))
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
                               SP_MtoD       MtoD)
  {
    SP_sweepsource p(new SweepSource(state, mesh, quadrature, material, MtoD));
    return p;
  }

  /*!
   * \brief Set an external moment source.
   *
   * \param source  Smart pointer to external source
   */
  void set_moment_source(SP_externalsource source)
  {
    Require(source);
    d_moment_external_sources.push_back(source);
  }

  /*!
   * \brief Set an external discrete source.
   *
   * \param source  Smart pointer to external source
   */
  void set_discrete_source(SP_externalsource source)
  {
    Require(source);
    d_discrete_external_sources.push_back(source);
  }

  /*!
   * \brief Set a fission source.
   *
   * \param source  Smart pointer to fission source
   */
  void set_fission_source(SP_fissionsource source)
  {
    Require(source);
    d_fissionsource = source;
  }

  // The following source construction routines aim to
  // satisfy the needs of all solvers.

  /*!
   *  \brief Add all fixed moments source.
   *
   *  Creates a fixed group source using external moments sources
   *  and, if applicable, fission sources.  Note, fission is
   *  only fixed within an eigenproblem.
   */
  void build_fixed(int g);

  /*!
   * \brief Build fixed with in-scatter.
   *
   * Creates a fixed group source as \ref build_fixed, but adds
   * the in-scatter source.  This is the typical construction
   * to be called within source iteration.
   *
   */
  void build_fixed_with_scatter(int g);

  /*!
   *  \brief Build fixed with downs-scatter.
   *
   *  Creates a fixed group source as \ref build_fixed, but adds
   *  the in-scatter source.  This is the typical construction
   *  to be called within source iteration.
   *
   *  \param g          Group of current solve
   *  \param g_cutoff   Last group to include in downscatter source
   */
  void build_fixed_with_downscatter(int g, int g_cutoff);

  /*!
   * \brief Build within-group scattering source.
   *
   * Called before each sweep.
   */
  void build_within_group_scatter(int g, const moments_type &phi);

  /*!
   *  \brief Build total scattering source.
   *
   *  This builds the complete group scatter source using the
   *  given multigroup flux vector.  This routine would find
   *  use in a multigroup Krylov solver.
   */
  void build_total_scatter(int g, int g_cutoff, const State::vec_moments_type &phi);

  /// Reset all the internal source vectors to zero.
  void reset();

  /*!
   *  \brief Get the sweep source vector.
   *
   */
  const sweep_source_type& source(int g, int o, int a);
  void source(int g, int o, int a, sweep_source_type& s);

private:

  SP_state d_state;
  SP_mesh d_mesh;
  SP_quadrature d_quadrature;
  SP_MtoD d_MtoD;

  /// Sweep source for a given angle and group over all cells.
  sweep_source_type d_source;

  /// Fixed moments source applicable to all angles in this group.
  moments_type d_fixed_group_source;

  /// Within group scattering applicable to all angles in this group.
  moments_type d_scatter_group_source;

  /// A container of external moments sources.
  std::vector<SP_externalsource> d_moment_external_sources;

  /// A container of external discrete sources
  std::vector<SP_externalsource> d_discrete_external_sources;

  /// Fission source
  SP_fissionsource d_fissionsource;

  /// Scattering source
  SP_scattersource d_scattersource;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "SweepSource.i.hh"

#endif /* SWEEPSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of SweepSource.hh
//---------------------------------------------------------------------------//
