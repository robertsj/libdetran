//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LRA.hh
 *  @author robertsj
 *  @date   Nov 29, 2012
 *  @brief  LRA class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_user_LRA_HH_
#define detran_user_LRA_HH_

#include "kinetics/TimeDependentMaterial.hh"
#include "geometry/Mesh.hh"

namespace detran_user
{

class LRA: public detran::TimeDependentMaterial
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran::TimeDependentMaterial             Base;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mesh       User-defined LRA mesh
   */
  LRA(SP_mesh mesh, bool flag, bool steady);

  // SP constructor
  static SP_material Create(SP_mesh, bool flag, bool steady);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  void set_state(SP_state);
  void initialize_materials();
  void update_P_and_T(double t);
  vec_dbl T() {return d_T;}
  vec_dbl P() {return d_P;}
  double area() {return d_A;}

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Mesh
  SP_mesh d_mesh;
  /// Diffusion flag
  bool d_flag;
  /// Steady state flag
  bool d_steady;
  /// Fine mesh map of unique coarse meshes (for base materials)
  vec_int d_unique_mesh_map;
  /// Fine mesh map of assemblies
  vec_int d_assembly_map;
  /// Fine mesh temperature
  vec_dbl d_T;
  /// Old fine mesh temperature
  vec_dbl d_T_old;
  /// Fine mesh power density
  vec_dbl d_P;
  /// Core area
  double d_A;
  /// Store the current time for iteration purposes
  double d_current_time;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL TIME DEPENDENT MATERIALS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /*
   *  @brief User-defined material update
   *
   *  This routine is called from within update.  This fills the
   *  internal cross section vectors with their actual values.  The
   *  update function then adds the synthetic components.
   *
   */
  void update_impl();

};

} // end namespace detran_user

#endif /* detran_user_LRA_HH_ */
