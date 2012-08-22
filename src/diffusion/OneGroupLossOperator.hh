//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OneGroupLossOperator.hh
 * \brief  OneGroupLossOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef ONEGROUPLOSSOPERATOR_HH_
#define ONEGROUPLOSSOPERATOR_HH_

// Detran
#include "BaseOperator.hh"

// Utilities
#include "Definitions.hh"

namespace detran_diffusion
{

/*!
 *  \class OneGroupLossOperator
 *  \brief Loss operator for a one group diffusion equation.
 *
 *  This operator can be used for group sweeping in multigroup
 *  Krylov diffusion solvers.  Additionally, it can be used
 *  to precondition within-group transport solves.
 *
 *  The one group diffusion equation is
 *  \f[
 *  - \nabla D(\vec{r}) \nabla \phi(\vec{r})
 *    + \Sigma_r(\vec{r}) \phi(\vec{r}) = Q(\vec{r})
 *  \f]
 *  where the loss operator is defined
 *  \f[
 *    \mathbf{M}[\cdot] \equiv (-\nabla D(\vec{r}) \nabla
 *      + \Sigma_r(\vec{r}))[\cdot]\, .
 *  \f]
 *
 *  A mesh-centered discretization is used.
 */
class OneGroupLossOperator: public BaseOperator
{

public:

  typedef detran::SP<OneGroupLossOperator>  SP_operator;
  typedef detran::InputDB::SP_input         SP_input;
  typedef detran::Material::SP_material     SP_material;
  typedef detran::Mesh::SP_mesh             SP_mesh;
  typedef detran::vec_dbl                   vec_dbl;

  /*!
   *  \brief Constructor
   *  \param
   */
  OneGroupLossOperator(SP_input    input,
                       SP_material material,
                       SP_mesh     mesh,
                       int         group);

private:

  /// \name Private Data
  /// \{

  /// Energy group for this operator
  int d_group;

  /// Albedos
  vec_dbl d_albedo;

  /// \}

  /// \name Implementation
  /// \{

  void construct();

  /// \}

};

} // end namespace detran

#endif // ONEGROUPLOSSOPERATOR_HH_ 

//---------------------------------------------------------------------------//
//              end of file OneGroupLossOperator.hh
//---------------------------------------------------------------------------//
