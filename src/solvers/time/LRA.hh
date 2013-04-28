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

// Raw material data    0                   1                   2                   3                   4                   5
//                      fuel 1, blade in    fuel 1, blade out   fuel 2, blade in    fuel 2, blade out   reflector           cr (= fuel 2, blade in @ t = 0)
const double T1[]  =  { 0.265604249667995 , 0.262881177707676 , 0.264760391845380 , 0.264760391845380 , 0.265181649429859 , 0.264760391845380 };
const double T2[]  =  { 1.579778830963665 , 1.752541184717841 , 1.594133588394707 , 1.594133588394707 , 2.093802345058626 , 1.594133588394707 };
const double A1[]  =  { 0.008252000000000 , 0.007181000000000 , 0.008002000000000 , 0.008002000000000 , 0.000602400000000 , 0.008002000000000 };
const double A2[]  =  { 0.100300000000000 , 0.070470000000000 , 0.083440000000000 , 0.073324000000000 , 0.019110000000000 , 0.083440000000000 };
const double F1[]  =  { 0.001893827160494 , 0.001896707818930 , 0.001918930041152 , 0.001918930041152 , 0.000000000000000 , 0.001918930041152 };
const double F2[]  =  { 0.044897119341564 , 0.035699588477366 , 0.042016460905350 , 0.042016460905350 , 0.000000000000000 , 0.042016460905350 };
const double S11[] =  { 0.232022249667995 , 0.228030177707676 , 0.230588391845380 , 0.230588391845380 , 0.21703924942985  , 0.230588391845380 };
const double S21[] =  { 0.025330000000000 , 0.027670000000000 , 0.026170000000000 , 0.026170000000000 , 0.047540000000000 , 0.026170000000000 };
const double S22[] =  { 1.479478830963665 , 1.682071184717841 , 1.510693588394707 , 1.520809588394707 , 2.074692345058626 , 1.510693588394707 };
const double D1[]  =  { 1.255000000000000 , 1.268000000000000 , 1.259000000000000 , 1.259000000000000 , 1.257000000000000 , 1.259000000000000 };
const double D2[]  =  { 0.211000000000000 , 0.190200000000000 , 0.209100000000000 , 0.209100000000000 , 0.159200000000000 , 0.209100000000000 };
const double mu0[] =  { 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 };
const double mu1[] =  { 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 , 0.000000000000000 };
const double NU = 2.43;
const double B  = 0.0001;
const double ALPHA   = 3.830e-11;
const double GAMMA   = 3.034e-3; // Benchmark spec has 2.034e-3
const double KAPPA   = 3.204e-11;
const double LAMBDA0 = 0.0654; // Benchmark spec has 0.00654
const double LAMBDA1 = 1.35;
const double BETA0   = 0.0054;
const double BETA1   = 0.0010873;
const double VELOCITY0 = 3.0e7;
const double VELOCITY1 = 3.0e5;
const int ROD = 5;
const int REFLECTOR = 4;

class SOLVERS_EXPORT LRA: public detran::TimeDependentMaterial
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
