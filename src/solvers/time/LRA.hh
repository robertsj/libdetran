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
#include "kinetics/MultiPhysics.hh"
#include "geometry/Mesh.hh"
#include "TimeStepper.hh"

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
  typedef detran::MultiPhysics::SP_multiphysics     SP_multiphysics;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mesh       User-defined LRA mesh
   */
  LRA(SP_mesh mesh, bool doingtransport, bool steady);

  // SP constructor
  static SP_material Create(SP_mesh, bool flag, bool steady);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  void set_state(SP_state);
  void initialize_materials();
  void update_P_and_T(double t, double dt);
  vec_dbl T() {return d_T;}
  vec_dbl P() {return d_P;}
  SP_multiphysics physics() {return d_physics;}
  double area() {return d_A;}
  void set_area(double a) {d_A = a;}

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
  /// Multiphysics (just T)
  SP_multiphysics d_physics;
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

//---------------------------------------------------------------------------//
template <class D>
void update_T_rhs(void* data,
                  detran::TimeStepper<D>* step,
                  double t,
                  double dt)
{
  Require(data);
  Require(step);

  // cast data as LRA
  LRA* mat = (LRA*) data;

  std::cout << mat->physics()->variable(0)[0] << " "
            << step->multiphysics()->variable(0)[0] << std::endl;
  // update
  mat->update_P_and_T(t, dt);
}
} // end namespace detran_user

#endif /* detran_user_LRA_HH_ */
