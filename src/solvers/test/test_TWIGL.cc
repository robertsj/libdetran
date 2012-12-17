//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_TWIGL.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  Implements the TWIGL kinetics benchmark
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST               \
        FUNC(test_TWIGL)

#include "TestDriver.hh"
#include "TimeStepper.hh"
#include "Mesh2D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "solver/EigenvalueManager.hh"
//
#include "angle/test/quadrature_fixture.hh"
#include "geometry/test/mesh_fixture.hh"
#include "material/test/material_fixture.hh"
#include "external_source/test/external_source_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_material;
using namespace detran_external_source;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

void test_monitor(void* data, TimeStepper<_1D>* ts, int step, double t,
                  double dt, int it, bool conv)
{
  printf(" %16.13f  %16.13f  %16.13f \n",
         t, ts->state()->phi(0)[0], ts->precursor()->C(0)[0]);
}

//---------------------------------------------------------------------------//
class TWIGLMaterial: public TimeDependentMaterial
{
public:
  typedef TimeDependentMaterial Base;
  TWIGLMaterial(bool reactivity = true)
    : Base(3, 2, 1, "TWIGLMaterial")
    , d_reactivity(reactivity)
  {
    /* ... */
  }
  // Update the materials.
  void update_impl()
  {
    double t = time();
    //
    int m = 0;
    set_diff_coef(m, 0,  1.400);
    set_diff_coef(m, 1,  0.400);

    double sa = 0.0;
    if (d_reactivity)
    {
      if (t < 0.2)
        sa = 0.150 * (1.0 - 0.11667 * t);
      else
        sa = 0.150 * 0.97666;
    }
    else
    {
      if (t == 0.0)
        sa = 0.150;
      else
        sa = 0.150 - 0.0;
    }
    set_sigma_t(m, 0,    0.010 + 0.010);
    set_sigma_t(m, 1,    sa);
    set_sigma_s(m, 1, 0, 0.010);
    set_sigma_f(m, 0,    0.007);
    set_sigma_f(m, 1,    0.200);
    set_chi(m, 0,        1.000);
    //
    m = 1;
    set_diff_coef(m, 0,  1.400);
    set_diff_coef(m, 1,  0.400);
    set_sigma_t(m, 0,    0.010 + 0.010);
    set_sigma_t(m, 1,    0.150);
    set_sigma_s(m, 1, 0, 0.010);
    set_sigma_f(m, 0,    0.007);
    set_sigma_f(m, 1,    0.200);
    set_chi(m, 0,        1.000);
    //
    m = 2;
    set_diff_coef(m, 0,  1.300);
    set_diff_coef(m, 1,  0.500);
    set_sigma_t(m, 0,    0.008 + 0.01);
    set_sigma_t(m, 1,    0.050);
    set_sigma_s(m, 1, 0, 0.010);
    set_sigma_f(m, 0,    0.003);
    set_sigma_f(m, 1,    0.060);
    set_chi(m, 0,        1.000);
    // kinetics
    set_lambda(0,        0.08);
    set_beta(0,          0.0*0.0075);
    set_velocity(0,      1.0e7);
    set_velocity(1,      2.0e5);
    // delayed chi(m, i, g, v)
    set_chi_d(0, 0, 0,    1.00000);
    set_chi_d(1, 0, 0,    1.00000);
    set_chi_d(2, 0, 0,    1.00000);
    //
    finalize();
  }

private:
  /// True for ramp; false for step.
  bool d_reactivity;
};

//---------------------------------------------------------------------------//
Mesh2D::SP_mesh get_mesh(Mesh2D::size_t fmm = 1)
{
  Mesh2D::vec_dbl cm(4);
  cm[0] =  0.0;
  cm[1] = 24.0;
  cm[2] = 56.0;
  cm[3] = 80.0;
  Mesh2D::vec_int fm(3);
  fm[0] = fmm * 3;
  fm[1] = fmm * 4;
  fm[2] = fmm * 3;
  Mesh2D::vec_int mt(9);
  mt[0] = 2; mt[1] = 1; mt[2] = 2;
  mt[3] = 1; mt[4] = 0; mt[5] = 2;
  mt[6] = 2; mt[7] = 2; mt[8] = 2;
  Mesh2D::SP_mesh mesh = Mesh2D::Create(fm, fm, cm, cm, mt);
  return mesh;
}

//---------------------------------------------------------------------------//
int test_TWIGL(int argc, char *argv[])
{

  typedef TimeStepper<_2D> TS_2D;

  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("TWIGL benchmark"));
  inp->put<int>("dimension",                2);
  inp->put<int>("number_groups",            2);
  inp->put<std::string>("equation",         "diffusion");
  inp->put<std::string>("bc_west",          "reflect");
  inp->put<std::string>("bc_east",          "vacuum");
  inp->put<std::string>("bc_south",         "reflect");
  inp->put<std::string>("bc_north",         "vacuum");
  inp->put<int>("bc_zero_flux",             1);
  inp->put<double>("ts_final_time",         1.0);
  inp->put<double>("ts_step_size",          0.1);
  inp->put<int>("ts_max_steps",             10);
  inp->put<int>("ts_scheme",                TS_2D::BDF1);
  inp->put<int>("ts_discrete",              0);
  inp->put<int>("ts_output",                0);
  inp->put<int>("ts_monitor_level",         1);
  // inner gmres parameters
  InputDB::SP_input db(new InputDB("inner_solver_db"));
  db->put<double>("linear_solver_atol",                 1e-12);
  db->put<double>("linear_solver_rtol",                 1e-12);
  db->put<string>("linear_solver_type",                 "gmres");
  db->put<int>("linear_solver_maxit",                   2000);
  db->put<int>("linear_solver_gmres_restart",           20);
  db->put<int>("linear_solver_monitor_level",           0);
  inp->put<InputDB::SP_input>("inner_solver_db",        db);
  inp->put<InputDB::SP_input>("outer_solver_db",        db);

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  // Create a TWIGL material with ramp reactivity
  TS_2D::SP_material mat(new TWIGLMaterial(true));

  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  TS_2D::SP_mesh mesh = get_mesh();

  //-------------------------------------------------------------------------//
  // STEADY STATE
  //-------------------------------------------------------------------------//

  EigenValueManager<_2D> eigen(inp, mat, mesh);
  FixedSourceManager<_1D> manager(inp, mat, mesh);
  manager.setup();
  manager.set_source(q_e1);
  manager.set_solver();
  //manager.solve();
  State::SP_state ic = manager.state();


  //-------------------------------------------------------------------------//
  // TIME STEPPER
  //-------------------------------------------------------------------------//


  TS_1D stepper(inp, linmat, mesh, true);
  stepper.set_monitor(test_monitor);

//  // Initial condition (constant psi = 1/2)
//  for (int o = 0; o < stepper.quadrature()->number_octants(); ++o)
//  {
//    for (int a = 0; a < stepper.quadrature()->number_angles_octant(); ++a)
//    {
//      for (int i = 0; i < mesh->number_cells(); ++i)
//      {
//        ic->phi(0)[i] = 1.0;
//        ic->psi(0, o, a)[i] = 0.5;
//      }
//    }
//  }
  // Initial condition (constant psi = 1/2)
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    ic->phi(0)[i] = 1.0;
  }

  //stepper.add_source(q_td);

  stepper.solve(ic);

  printf(" %20.16f %20.16f ", ic->phi(0)[0], ic->phi(0)[1]);
  std::cout << std::endl;

  State::SP_state final = stepper.state();

  printf(" %20.16f %20.16f ", final->phi(0)[0], final->phi(0)[1]);
  std::cout << std::endl;

  return 0;

}
//---------------------------------------------------------------------------//
//              end of test_TWIGL.cc
//---------------------------------------------------------------------------//
