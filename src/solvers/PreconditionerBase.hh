//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerBase.hh
 * \brief  PreconditionerBase 
 * \author Jeremy Roberts
 * \date   Aug 5, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONERBASE_HH_
#define PRECONDITIONERBASE_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "ScatterSource.hh"

// Diffusion
#include "BaseOperator.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"

// System
#include "petsc.h"

namespace detran
{

class PreconditionerBase
{

public:

  typedef InputDB::SP_input               SP_input;
  typedef Material::SP_material           SP_material;
  typedef Mesh::SP_mesh                   SP_mesh;
  typedef ScatterSource::SP_source        SP_scattersource;
  typedef detran_diffusion::
          BaseOperator::SP_operator       SP_lossoperator;

  /*!
   *  \class PreconditionerBase
   *  \brief Base class for defining a preconditioner.
   *
   * Both the within-group and multigroup transport problems
   * can be put in the form
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
   * Both a within-group (\ref PreconditionerWG) and a multigroup
   * (\ref PreconditionerMG) are implemented using the
   * associated within-group diffusion operator
   * (\ref OneGroupLossOperator) and multigroup diffusion
   * operator (\ref LossOperator).
   *
   * In both cases, the preconditioning process \$\mathbf{P}^{-1} \$
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
  PreconditionerBase(SP_input input,
                     SP_material material,
                     SP_mesh mesh,
                     SP_scattersource scattersource)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_scattersource(scattersource)
  {
    Require(d_input);
    Require(d_material);
    Require(d_mesh);
    Require(d_scattersource);
  }

  /// Virtual destructor
 virtual ~PreconditionerBase(){};

 /*!
  *  \brief Apply the preconditioning process, \f$ \mathbf{P}^{-1} \f$.
  */
 virtual PetscErrorCode apply(Mat A, Vec x, Vec y) = 0;

protected:

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
  int d_size;

  /// Preconditioner
  Mat d_operator;

  /// Linear solver
  KSP d_solver;

  /// Diffusion loss operator
  SP_lossoperator d_lossoperator;

  /// \}

};

} // end namespace detran

#endif // PRECONDITIONERBASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file PreconditionerBase.hh
//---------------------------------------------------------------------------//
