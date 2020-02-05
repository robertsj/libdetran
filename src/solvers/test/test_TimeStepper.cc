//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_TimeStepper.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  test_TimeStepper class definition.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_TimeStepper)      \
        FUNC(test_BDF_Steps)

#include "TestDriver.hh"
#include "TimeStepper.hh"
#include "Mesh1D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
#include "kinetics/LinearMaterial.hh"
#include "external_source/IsotropicSource.hh"

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

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

void test_monitor(void* data, TimeStepper<_1D>* ts, int step, double t, double dt, int it, bool conv)
{
  printf(" %16.13f  %16.13f  %16.13f \n", t, ts->state()->phi(0)[0], ts->precursor()->C(0)[0]);
}

/**
 *  1D time dependent SN problem with a TD source and no
 *  fission.  Constant source from 0:1.  Final time 10.
 */
int test_TimeStepper(int argc, char *argv[])
{

  typedef TimeStepper<_1D> TS_1D;

  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("time stepper test"));
  inp->put<int>("dimension",                1);
  inp->put<int>("number_groups",            1);
  inp->put<std::string>("equation",         "diffusion");
  inp->put<std::string>("bc_west",          "reflect");
  inp->put<std::string>("bc_east",          "reflect");
  inp->put<double>("ts_final_time",         1.0); // 0.95162581964040427
  inp->put<double>("ts_step_size",          0.01);
  inp->put<int>("ts_max_steps",             1000);
  inp->put<int>("ts_scheme",                TS_1D::IMP);
  inp->put<int>("ts_discrete",              0);
  inp->put<int>("ts_output",                0);
  inp->put<int>("ts_monitor_level",         1);
  inp->put<int>("store_angular_flux",       0);
  inp->put<int>("inner_print_level",        0);
  inp->put<int>("outer_print_level",        0);
  inp->put<int>("quad_number_polar_octant", 1);
  inp->put<double>("inner_tolerance",       1e-19);
  inp->put<int>("inner_max_iters",          1e7);
  inp->put<int>("compute_boundary_flux",    1);
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

  // One group, homogeneous, and time-invariant.
  // We need to create the KineticsMaterial, and use that
  // as the one material in a LinearMaterial.

  KineticsMaterial::SP_material
    kinmat(new KineticsMaterial(1, 1, 1, "base material"));
  kinmat->set_sigma_t(0, 0,    1.0);
  kinmat->set_diff_coef(0, 0,  1.0/3.0);
  //kinmat->set_sigma_s(0, 0, 0, 0.5);
  kinmat->set_sigma_f(0, 0,    0.0);
  kinmat->set_velocity(0,      1.0);
  kinmat->set_chi(0, 0,        1.0);
  kinmat->set_beta(0, 0,       0.00);
  kinmat->set_chi_d(0, 0, 0,   1.0);
  kinmat->set_lambda(0,        0.1);
  kinmat->finalize();
  LinearMaterial::vec_material materials(1, kinmat);
  LinearMaterial::vec_dbl      times(1, 0.0);
  LinearMaterial::SP_material  linmat(new LinearMaterial(times, materials));

  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  vec_int fm(1, 3);
  vec_dbl cm(2, 0.0); cm[1] = 1.0;
  vec_int mt(1, 0);
  Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));

  //-------------------------------------------------------------------------//
  // SOURCE
  //-------------------------------------------------------------------------//

  // We want a source in the left half of
  // the slab, turned on from 0 to 1 seconds.  We half
  // to approximate via a small ramp between 1 and 1 + eps.
  //
  // Static source
  vec_int source_map(mesh->number_cells(), 0);
  for (int i = 0; i < mesh->number_cells(); ++i)
    source_map[i] = 1;
  IsotropicSource::spectra_type spectra(2, vec_dbl(1, 0.0));
  // "on"
  spectra[1][0] = 1.0;
  IsotropicSource::SP_externalsource
    q_e1(new IsotropicSource(1, mesh, spectra, source_map));
  // "off"
  spectra[1][0] = 1.0;
  IsotropicSource::SP_externalsource
    q_e2(new IsotropicSource(1, mesh, spectra, source_map));
  //
  // Linear source
  double time_off = 100.0;
  vec_dbl source_times(2, time_off); source_times[1] = time_off + 0.0;
  LinearExternalSource::vec_source sources;
  sources.push_back(q_e1);
  sources.push_back(q_e2);
  LinearExternalSource::SP_tdsource
    q_td(new LinearExternalSource(1, mesh, source_times, sources));

  //-------------------------------------------------------------------------//
  // INITIAL CONDITION
  //-------------------------------------------------------------------------//

  // Here, we want the flux due to the source being "on" for all time
  // before t = 0.  This means we need to solve the fixed source problem
  // using that initial source.

  FixedSourceManager<_1D> manager(inp, kinmat, mesh);
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

  // Initial condition (constant psi = 1/2)
  for (int i = 0; i < mesh->number_cells(); ++i)
  {
    ic->phi(0)[i] = 1.0;
  }

  stepper.add_source(q_td);

  stepper.solve(ic);

  printf(" %20.16f %20.16f ", ic->phi(0)[0], ic->phi(0)[1]);
  std::cout << std::endl;

  State::SP_state final = stepper.state();

  printf(" %20.16f %20.16f ", final->phi(0)[0], final->phi(0)[1]);
  std::cout << std::endl;

  return 0;

}

/// Test BDF steps agains reference
int test_BDF_Steps(int argc, char *argv[])
{

  typedef TimeStepper<_1D> TS_1D;

  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("time stepper test"));
  inp->put<int>("dimension",                1);
  inp->put<int>("number_groups",            1);
  inp->put<std::string>("equation",         "sc");
  inp->put<std::string>("bc_west",          "reflect");
  inp->put<std::string>("bc_east",          "reflect");
  inp->put<double>("ts_final_time",         1.0);
  inp->put<double>("ts_step_size",          0.1);
  inp->put<int>("ts_discrete",              1);
  inp->put<int>("ts_output",                0);
  inp->put<int>("ts_monitor_level",         0);
  inp->put<int>("store_angular_flux",       1);
  inp->put<int>("inner_print_level",        0);
  inp->put<int>("outer_print_level",        0);
  inp->put<int>("quad_number_polar_octant", 1);
  inp->put<double>("inner_tolerance",       1e-16);
  inp->put<int>("inner_max_iters",          1e7);
  inp->put<int>("compute_boundary_flux",    1);

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  // Time-invariant cross-section with sigmaT=sigmaA=1.0
  KineticsMaterial::SP_material
    kinmat(new KineticsMaterial(1, 1, 0, "base material"));
  kinmat->set_sigma_t(0, 0,    1.0);
  kinmat->set_velocity(0,      1.0);
  kinmat->finalize();
  LinearMaterial::vec_material materials(1, kinmat);
  LinearMaterial::vec_dbl      times(1, 0.0);
  LinearMaterial::SP_material  linmat(new LinearMaterial(times, materials));

  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  vec_int fm(1, 3);
  vec_dbl cm(2, 0.0); cm[1] = 1.0;
  vec_int mt(1, 0);
  Mesh1D::SP_mesh mesh(new Mesh1D(fm, cm, mt));

  //-------------------------------------------------------------------------//
  // TEST BDF
  //-------------------------------------------------------------------------//

  // Reference fluxes
  double ref[] = {0.9047619047619050,  // IMP   @ 0.1
                  0.9090909090909091,  // BDF1  @ 0.1
                  0.8184523809523810,  // BDF2  @ 0.2
                  0.7404556650246303,  // BDF3  @ 0.3
                  0.6699914764537027,  // BDF4  @ 0.4
                  0.6062522850493363,  // BDF5  @ 0.5
                  0.5485489911896624}; // BDF6  @ 0.6


  for (int scheme = 0; scheme < 7; ++scheme)
  {
    int max_steps = scheme;
    if (max_steps == 0) max_steps = 1;
    inp->put<int>("ts_max_steps",             max_steps);
    inp->put<int>("ts_scheme",                scheme);

    // Stepper
    TS_1D stepper(inp, linmat, mesh, false);

    // Initial condition (constant psi = 1/2)
    TS_1D::SP_state ic = stepper.state();
    for (int o = 0; o < stepper.quadrature()->number_octants(); ++o)
      for (int a = 0; a < stepper.quadrature()->number_angles_octant(); ++a)
        for (int i = 0; i < mesh->number_cells(); ++i)
          ic->psi(0, o, a)[i] = 0.5;

    stepper.solve(ic);

    State::SP_state final = stepper.state();
    std::cout << " phi = " << final->phi(0)[0]
              << " ref = " << ref[scheme] << std::endl;
    TEST(soft_equiv(final->phi(0)[0], ref[scheme]));

  } // end scheme loop

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TimeStepper.cc
//---------------------------------------------------------------------------//
