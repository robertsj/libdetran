//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGScatterFissionOperator.hh
 *  @brief MGScatterFissionOperator class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_MGSCATTERFISSIONOPERATOR_HH_
#define detran_MGSCATTERFISSIONOPERATOR_HH_

#include "callow/matrix/MatrixShell.hh"
#include "transport/FissionSource.hh"
#include "transport/ScatterSource.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"

namespace detran
{

/**
 *  @class MGScatterFissionOperator
 *  @brief Wraps the scatter (and fission) source in operator form
 */
class MGScatterFissionOperator: public callow::MatrixShell
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum SFOPTIONS
  {
    SCATTERONLY, SCATTERFISSION, FISSIONONLY, END_SFOPTIONS
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MatrixShell                                 Base;
  typedef MGScatterFissionOperator                    Operator_T;
  typedef detran_utilities::SP<Operator_T>            SP_operator;
  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef detran_material::Material::SP_material      SP_material;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef ScatterSource::SP_scattersource             SP_scattersource;
  typedef FissionSource::SP_fissionsource             SP_fissionsource;
  typedef detran_utilities::size_t                    size_t;
  typedef detran_utilities::vec_size_t                groups_t;
  typedef groups_t::iterator                          groups_iter;
  typedef State::moments_type                         moments_type;
  typedef callow::Vector                              Vector;

  /**
   *  @brief Constructor
   *  @param    input     parameter database
   *  @param    material  cross sections
   *  @param    mesh      geometry
   *  @param    S         scatter source
   *  @param    F         fission source
   *  @param    cutoff    first group included in solve
   *  @param    sf_switch 0 for scatter only, 1 for both, 2 for fission only
   *  @param    adjoint   flag for adjoint problems
   */
  MGScatterFissionOperator(SP_input         input,
                           SP_material      material,
                           SP_mesh          mesh,
                           SP_scattersource S,
                           SP_fissionsource F,
                           size_t           cutoff,
                           size_t           sf_switch,
                           bool             adjoint);

  /// Virtual destructor
  virtual ~MGScatterFissionOperator(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  // the client must implement the action y <-- A * x
  void multiply(const Vector &x,  Vector &y);

  // the client must implement the action y <-- A' * x
  void multiply_transpose(const Vector &x, Vector &y)
  {
    THROW("Transpose scatter+fission transpose operator not implemented");
  }
private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_input d_input;
  SP_material d_material;
  SP_mesh d_mesh;
  SP_scattersource d_S;
  SP_fissionsource d_F;
  /// Cutoff
  size_t d_group_cutoff;
  /// Scatter flag
  bool d_include_scatter;
  /// Fission flag
  bool d_include_fission;
  /// Adjoint flag
  bool d_adjoint;
  /// Number of groups
  size_t d_number_groups;
  /// Number of active groups
  size_t d_number_active_groups;
  size_t d_moments_size;
  /// Groups to sweep
  groups_t d_groups;

};

} // end namespace detran

#endif /* detran_MGSCATTERFISSIONOPERATOR_HH_ */

//----------------------------------------------------------------------------//
//              end of file MGScatterFissionOperator.hh
//----------------------------------------------------------------------------//
