//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerMG.hh
 * \brief  PreconditionerMG
 * \author Jeremy Roberts
 * \date   Aug 6, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONERMG_HH_
#define PRECONDITIONERMG_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "ScatterSource.hh"

// Diffusion
#include "LossOperator.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"
#include "SP.hh"

namespace detran
{

/*!
 *  \class PreconditionerMG
 *  \brief Multigroup preconditioner.
 *
 * The multigroup transport problem can be put in the form
 * \f[
 *     \overbrace{(\mathbf{I} -
 *                \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})}^{\mathbf{A}} \phi
 *     = \overbrace{\mathbf{D} \mathbf{L}^{-1} Q}^{b} \, ,
 * \f]
 *
 * A good preconditioner \f$ \mathbf{P} \f$ for this problem
 * should satisfy
 * \f[
 *     \mathbf{P}^{-1} \mathbf{A} \approx \mathbf{I} \, ,
 * \f]
 * and, importantly, the action of \$\mathbf{P}^{-1} \$ should be
 * cheap to compute.
 *
 * Currently, we provide diffusion-based preconditioner that includes
 * all in-scatter.
 *
 * The preconditioning process \$\mathbf{P}^{-1} \$
 * is defined to be
 * \f[
 *     (\mathbf{I} - \mathbf{C}^{-1} \mathbf{S}) \, ,
 * \f]
 * where \f$ \mathbf{C} \f$ is the diffusion operator.
 *
 * \note The diffusion operator only operates on the scalar
 *       flux, and so addition of higher order moments will
 *       require restriction and projection operations.
 *
 */
class PreconditionerMG
{

public:

  typedef SP<PreconditionerMG>            SP_pc;
  typedef InputDB::SP_input               SP_input;
  typedef Material::SP_material           SP_material;
  typedef Mesh::SP_mesh                   SP_mesh;
  typedef ScatterSource::SP_source        SP_scattersource;
  typedef detran_diffusion::
    LossOperator::SP_operator             SP_lossoperator;

  /*!
   *  \brief Constructor
   *
   *  Assuming the within-group transport problem is set up,
   *  a KSP object exists from which the PC is extracted.  This
   *  PC is passed here to be constructed and for its application
   *  operator to be assigned.
   *
   *  \param input                Input database
   *  \param material             Material database
   *  \param mesh                 Cartesian mesh
   *  \param source               Scattering source
   *  \param moments_size_group   Unknown size for group moment vector
   *  \param boundary_size_group  Unknown size for group boundary vector
   *  \param mg_pc                PETSc PC object to be constructed
   */
  PreconditionerMG(SP_input input,
                   SP_material material,
                   SP_mesh mesh,
                   SP_scattersource source,
                   int moments_size_group,
                   int boundary_size_group,
                   int upscatter_cutoff,
                   PC mg_pc);

  /// Destructor
  ~PreconditionerMG()
  {
    KSPDestroy(&d_solver);
    VecDestroy(&d_x);
    VecDestroy(&d_y);
  }

  /*!
   *  \brief Apply the preconditioning process, \f$ \mathbf{P}^{-1} \f$.
   */
  PetscErrorCode apply(Vec x, Vec y);

private:

  /// \name Private Data
  /// \{

  /// Input database
  SP_input d_input;

  /// Material database
  SP_material d_material;

  /// Cartesian mesh
  SP_mesh d_mesh;

  /// Scatter source
  SP_scattersource d_scattersource;

  /// System size
  u_int d_size;

  /// Vector of linear solvers for applying the inverse diffusion operator
  KSP d_solver;

  /// Vector of diffusion loss operators for each group
  SP_lossoperator d_lossoperator;

  /// Vec of solver size for use in diffusion solve
  Vec d_x;
  Vec d_y;

  /// Size of the moments in a group
  u_int d_moments_size_group;

  /// Size of the boundary in a group
  u_int d_boundary_size_group;

  /// Upscatter cutoff (groups below and including don't get upscatter)
  u_int d_upscatter_cutoff;

  /// Tolerance for the preconditioner solve
  double d_tolerance;

  /// \}

};

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

/// Apply the PC
PetscErrorCode apply_inv_P_MG(PC pc, Vec x, Vec y);

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "PreconditionerMG.i.hh"

#endif // PRECONDITIONERMG_HH_

//---------------------------------------------------------------------------//
//              end of file PreconditionerMG.hh
//---------------------------------------------------------------------------//
