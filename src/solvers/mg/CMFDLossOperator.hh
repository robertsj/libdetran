//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CMFDLossOperator.hh
 *  @brief CMFDLossOperator class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_CMFDLossOperator_HH_
#define detran_CMFDLossOperator_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include "callow/matrix/Matrix.hh"
#include "transport/CurrentTally.hh"

namespace detran
{

/**
 *  @class CMFDLossOperator
 *  @brief Loss operator for multigroup CMFD acceleration
 */
template <class D>
class CMFDLossOperator: public callow::Matrix
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef callow::Matrix                                Base;
  typedef callow::MatrixBase::SP_matrix                 SP_matrix;
  typedef detran_utilities::SP<CMFDLossOperator>        SP_lossoperator;
  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_material::Material::SP_material        SP_material;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_geometry::Mesh                         Mesh;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;
  typedef detran_utilities::vec3_dbl                    vec3_dbl;
  typedef detran_utilities::vec_size_t                  groups_t;
  typedef groups_t::iterator                            groups_iter;
  typedef typename CurrentTally<D>::SP_currenttally     SP_tally;

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
  CMFDLossOperator(SP_input      input,
                   SP_material   material,
                   SP_mesh       mesh,
                   SP_tally      tally,
                   const bool    include_fission,
                   const bool    adjoint = false,
                   const double  keff = 1.0);

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
   */
  void construct(const vec2_dbl &phi,
                 double          keff = 1.0,
                 SP_material     mat  = SP_material(0),
                 bool            init = true);

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
  /// Current tally
  SP_tally d_tally;
  /// Problem dimension
  size_t d_dimension;
  /// Energy-dependent albedos
  vec2_dbl d_albedo;
  /// Energy groups
  size_t d_number_groups;
  /// Group indices
  groups_t d_groups;
  /// One group spatial size
  size_t d_group_size;
  /// Including fission?
  bool d_include_fission;
  /// Scaling factor for fission source
  double d_keff;
  /// Adjoint flag
  bool d_adjoint;
  /// Correct coupling coefficient to keep positive definiteness
  bool d_correct;
  /// Coupling coefficients
  vec3_dbl d_d_hat;
  double d_alpha;
  bool d_initial;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  // All matrix operators need this.  Here, we have the client call construct.
  void build(const vec2_dbl &phi);

  // Coupling coefficient
  //double couple

};

} // end namespace detran

#endif // detran_CMFDLossOperator_HH_

//----------------------------------------------------------------------------//
//              end of file CMFDLossOperator.hh
//----------------------------------------------------------------------------//
