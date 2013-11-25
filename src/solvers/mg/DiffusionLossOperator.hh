//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  DiffusionLossOperator.hh
 *  @brief DiffusionLossOperator
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_DIFFUSIONLOSSOPERATOR_HH_
#define detran_DIFFUSIONLOSSOPERATOR_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include "callow/matrix/Matrix.hh"

namespace detran
{

/**
 *  @class DiffusionLossOperator
 *  @brief Loss operator for multigroup diffusion problems
 *
 *  The loss operator represents neutron losses through interactions
 *  and leakage from a volume.  In Detran, we employ a mesh-centered
 *  finite difference approximation.  Moreover, we use an energy
 *  block structure so that the diffusion operators match with the
 *  underlying transport structure, which is useful for preconditioning
 *  applications.
 *
 *  The loss operator can be constructed with or without a fission
 *  contribution.  For multiplying fixed source problems, it is
 *  numerically more efficient to account for fission implicitly
 *  by bringing it to the left hand side.  However, there are cases
 *  where performing fission iteration is warranted, in which case
 *  not including the fission source is required.
 *
 */
class DiffusionLossOperator: public callow::Matrix
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef callow::Matrix                                Base;
  typedef callow::MatrixBase::SP_matrix                 SP_matrix;
  typedef detran_utilities::SP<DiffusionLossOperator>   SP_lossoperator;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_geometry::Mesh                         Mesh;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;
  typedef detran_utilities::vec_size_t                  groups_t;
  typedef groups_t::iterator                            groups_iter;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input            Pointer to input parameters
   *  @param material         Pointer to materials
   *  @param mesh             Pointer to mesh
   *  @param include_fission  Flag for including fission source implicitly
   *  @param cutoff           Lowest group to include in the operator
   *  @param adjoint          Adjoint flag
   *  @param keff             Fission scaling factor
   */
  DiffusionLossOperator(SP_input      input,
                        SP_material   material,
                        SP_mesh       mesh,
                        const bool    include_fission,
                        const size_t  cutoff = 0,
                        const bool    adjoint = false,
                        const double  keff = 1.0,
                        const size_t  sf_flag = 0);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Rebuild the matrix for a new fission scaling constant
   *
   *  This allows the client to rebuild the matrix after initial
   *  construction.  This is useful for response function generation
   *  as a function of keff or for time-dependent problems in which
   *  the pseudo-coefficients changes with time.
   *
   *  @param keff   Scaling parameter for fission source
   *  @param flag   Optionally omit any scatter and fission terms
   */
  void construct(double keff = 1.0, const size_t flag = 0);

  double albedo(const size_t side, const size_t g) const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Input
  SP_input d_input;
  /// Material
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// Problem dimension
  size_t d_dimension;
  /// Energy-dependent albedos
  vec2_dbl d_albedo;
  /// Energy groups
  size_t d_number_groups;
  /// Cutoff
  size_t d_group_cutoff;
  /// Group indices
  groups_t d_groups;
  /// Active groups
  size_t d_number_active_groups;
  /// One group spatial size
  size_t d_group_size;
  /// Including fission?
  bool d_include_fission;
  /// Scaling factor for fission source
  double d_keff;
  /// Adjoint flag
  bool d_adjoint;
  /// Skip scatter and fission
  size_t d_sf_flag;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  // All matrix operators need this.  Here, we have the client call construct.
  void build();

};

} // end namespace detran

#endif // detran_DIFFUSIONLOSSOPERATOR_HH_

//----------------------------------------------------------------------------//
//              end of file DiffusionLossOperator.hh
//----------------------------------------------------------------------------//
