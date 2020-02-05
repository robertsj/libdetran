//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCoarseMeshPreconditioner.hh
 *  @brief MGCoarseMeshPreconditioner class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGCOARSEMESHPRECONDITIONER_HH_
#define detran_MGCOARSEMESHPRECONDITIONER_HH_

#include "MGPreconditioner.hh"
#include "transport/CoarseMesh.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/matrix/Matrix.hh"

namespace detran
{

/**
 *  @class MGCoarseMeshPreconditioner
 *  @brief A base preconditioner for processes using a coarse space/energy mesh
 *
 *  When condensing to a coarse mesh, the cross section data needs to be
 *  averaged over space and energy in some way.  In all coarse mesh
 *  preconditioning schemes, the preconditioner operates on a fine mesh
 *  vector @f$ u_h \in \mathcal{R}^{n} @f$, where @f$ n @f$ is the fine
 *  mesh dimension.  In general, this action has
 *  the form
 *  @f[
 *    u_h \leftarrow P M R u_h \, ,
 *
 *  @f]
 *  where
 *  @f$R \in  \mathbb{R}^{m\times n} @f$ is the restriction operator,
 *  @f$P \in  \mathbb{R}^{n\times m} @f$ is the prolongation operator,
 *  @f$M \in  \mathbb{R}^{m\times m} @f$ is the coarse mesh preconditioning
 *  operator, and @f$ m @f$ is the coarse mesh dimension.  Derived
 *  classes provide implementations of @f$ M @f$.
 *
 *  @section mgcmpc_space Spatial Condensation
 *
 *  Restriction in space is based on a user-specified level defining the
 *  number of fine cells per coarse cell in each dimension. Averaging over
 *  space is based on volume-integration with an assumed space-dependent
 *  spectrum.
 *
 *
 *  @section mgcmpc_energy Energy Condensation
 *
 *  Restriction in energy is based on a user-specified coarse group map.
 *  An energy-dependent spectrum is required.
 *
 *  The most straightforward approach is to assume a unity flux for all
 *  groups, which leads to pure volume weighting in space and averaging
 *  in energy.
 *
 *  A second approach is to use a representative flux for condensation. The
 *  method used here is similar to the two-grid transport acceleration of
 *  Adams and Morel.  It employs a material-specific spectrum corresponding
 *  to the spectral radius of  Gauss-Seidel iteration in an infinite medium.
 *  Mathematically, this means solving
 *  @f[
 *    (\mathbf{T}-\mathbf{S}_L-\mathbf{S}_D)^{-1}\mathbf{S}_U \varphi
 *      = \rho \varphi \, ,
 *  @f]
 *  where @f$ \mathbf{T} @f$ contains the total cross sections, and
 *  @f$ \mathbf{S}_L @f$, @f$ \mathbf{S}_U @f$, and @f$ \mathbf{S}_D @f$
 *  contain the lower triangle, upper triangle, and diagonal components of
 *  the scattering matrix, respectively.  Note, when the problem is
 *  multiplying, the @f$ \mathbf{S} @f$ components also contain the
 *  fission transfer matrix.
 */

class MGCoarseMeshPreconditioner: public MGPreconditioner
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum condensation_options
  {
    CONDENSE_WITH_STATE,            // fine mesh/group state
    CONDENSE_WITH_FLAT_SPECTRUM,    // flat flux
    CONDENSE_WITH_GS_SPECTRUM,      // material-dependent gauss-seidel mode
    CONDENSE_WITH_PINCELL_SPECTRUM, // region-dependent pincell spectrum
    CONDENSE_WITH_FS_SPECTRUM,      // material-dependent flat fixed source mode
    CONDENSE_WITH_USER_SPECTRUM,    // user-defined spectrum via db
    END_CONDENSATION_OPTIONS
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef FissionSource::SP_fissionsource   SP_fissionsource;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef CoarseMesh::SP_coarsemesh         SP_coarsemesh;
  typedef callow::Matrix::SP_matrix         SP_matrix;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec2_dbl        vec2_dbl;
  typedef detran_utilities::vec_size_t      groups_t;
  typedef groups_t::iterator                groups_iter;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input            Input database
   *  @param material         Material database
   *  @param mesh             Cartesian mesh
   *  @param cutoff           First group included in solve
   *  @param include_fission  Flag to reat fission like scatter
   */
  MGCoarseMeshPreconditioner(SP_input         input,
                             SP_material      material,
                             SP_mesh          mesh,
                             SP_scattersource ssource,
                             SP_fissionsource fsource,
                             size_t           cutoff,
                             bool             include_fission,
                             bool             adjoint,
                             std::string      name = "MG-PC");

  /// Virtual destructor
  virtual ~MGCoarseMeshPreconditioner();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE --- ALL MULTIGROUP PRECONDITIONERS MUST IMPLEMENT
  //--------------------------------------------------------------------------//

  /// Solve Px = b
  virtual void apply(Vector &b, Vector &x) = 0;

  /**
   *  @brief Setup base components of the preconditioner
   *
   *  Using a spectrum or state vector, this base implementation condenses
   *  the materials on a coarse mesh in space and energy, producing the
   *  restriction and prolongation operators.  The derived implementation
   *  should call this if the default coarse meshing is desired.
   */
  virtual void build(const double keff = 1.0, SP_state state = SP_state(0)) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Material condensation option
  size_t d_condensation_option;
  /// Active groups
  groups_t d_groups;
  /// Fine group counts
  vec_int d_fine_per_coarse;
  /// Fine to coarse group map
  groups_t d_f2c_group_map;
  /// State for condensation
  SP_state d_state;
  /// Spectrum for condensation
  vec2_dbl d_spectrum;
  /// Spectrum for condensation
  vec2_dbl d_coarse_spectrum;
  /// Spectrum coarse mesh key
  std::string d_key;
  /// Coarse mesh
  SP_mesh d_coarsemesh;
  /// Coarse material
  SP_material d_c_material;
  /// Fine mesh active problem size
  size_t d_size_fine;
  /// Coarse mesh active problem size
  size_t d_size_coarse;
  /// Restriction operator
  SP_matrix d_restrict;
  /// Prolongation operator
  SP_matrix d_prolong;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  void build_restrict();
  void build_prolong();

  /**
   *  @brief Get the spectral shape factor in a fine mesh/group
   *  @param    i_f     fine mesh index
   *  @param    g_f     fine group index
   *  @param    i_c     coarse mesh index
   *  @param    g_c     coarse group index
   */
  double shape(const size_t i_f, const size_t g_f,
               const size_t i_c, const size_t g_c);
};

} // end namespace detran

#endif /* detran_MGCOARSEMESHPRECONDITIONER_HH_ */

//----------------------------------------------------------------------------//
//              end of file MGCoarseMeshPreconditioner.hh
//----------------------------------------------------------------------------//
