//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerWG.hh
 * \brief  PreconditionerWG
 * \author Jeremy Roberts
 * \date   Aug 5, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONERWG_HH_
#define PRECONDITIONERWG_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "ScatterSource.hh"

// Diffusion
#include "OneGroupLossOperator.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"
#include "SP.hh"

namespace detran
{

/*!
 *  \class PreconditionerWG
 *  \brief Withing-group preconditioner.
 *
 * The within-group transport problem can be put in the form
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
 * Currently, we provide diffusion-based preconditioning that is
 * in essence equivalent to diffusion synthetic acceleration.
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
class PreconditionerWG
{

public:

  typedef SP<PreconditionerWG>            SP_pc;
  typedef InputDB::SP_input               SP_input;
  typedef Material::SP_material           SP_material;
  typedef Mesh::SP_mesh                   SP_mesh;
  typedef ScatterSource::SP_source        SP_scattersource;
  typedef detran_diffusion::
    OneGroupLossOperator::SP_operator     SP_lossoperator;
  typedef std::vector<KSP>                vec_KSP;
  typedef std::vector<SP_lossoperator>    vec_lossoperator;


  /*!
   *  \brief Constructor
   *
   *  Assuming the within-group transport problem is set up,
   *  a KSP object exists from which the PC is extracted.  This
   *  PC is passed here to be constructed and for its application
   *  operator to be assigned.
   *
   *  \param input      Input database
   *  \param material   Material database
   *  \param mesh       Cartesian mesh
   *  \param source     Scattering source
   *  \param wg_pc      PETSc PC object to be constructed
   */
  PreconditionerWG(SP_input input,
                   SP_material material,
                   SP_mesh mesh,
                   SP_scattersource source,
                   PC wg_pc);

  /// Destructor
  ~PreconditionerWG()
  {
    for (int g = 0; g < d_material->number_groups(); g++)
    {
      KSPDestroy(&d_solver[g]);
    }
    VecDestroy(&d_x);
    VecDestroy(&d_y);
  }

  /// Set the group for this solve.
  void set_group(u_int group)
  {
    // Preconditions
    Require(group < d_material->number_groups());
    d_group = group;
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
  vec_KSP d_solver;

  /// Vector of diffusion loss operators for each group
  vec_lossoperator d_lossoperator;

  /// Vec of solver size for use in diffusion solve
  Vec d_x;
  Vec d_y;

  /// Group of within group solve
  u_int d_group;

  /// Tolerance for the preconditioner solve
  double d_tolerance;

  /// \}

};

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

/// Apply the PC
PetscErrorCode apply_inv_P(PC pc, Vec x, Vec y);

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "PreconditionerWG.i.hh"

#endif // PRECONDITIONERWG_HH_

//---------------------------------------------------------------------------//
//              end of file PreconditionerWG.hh
//---------------------------------------------------------------------------//
