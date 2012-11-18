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
        FUNC(test_TimeStepper)

#include "TestDriver.hh"
#include "TimeStepper.hh"
#include "Mesh1D.hh"
#include "external_source/ConstantSource.hh"
#include "callow/utils/Initialization.hh"
#include "callow/vector/Vector.hh"
#include "boundary/BoundaryTraits.hh"
#include "boundary/BoundarySN.hh"
#include "kinetics/LinearExternalSource.hh"
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

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

/**
 *  1D time dependent SN problem with a TD source and no
 *  fission.  Constant source from 0:1.  Final time 10.
 */
int test_TimeStepper()
{
  //-------------------------------------------------------------------------//
  // INPUT
  //-------------------------------------------------------------------------//

  InputDB::SP_input inp(new InputDB("time stepper test"));
  inp->put<int>("dimension",          1);
  inp->put<int>("number_groups",      1);
  inp->put<double>("ts_final_time",   10.0);
  inp->put<double>("ts_step_size",    1.0);
  inp->put<int>("ts_scheme",          BDF1);
  inp->put<int>("ts_discrete",        1);

  //-------------------------------------------------------------------------//
  // MATERIAL
  //-------------------------------------------------------------------------//

  // One group, homogeneous, and time-invariant.
  // We need to create the KineticsMaterial, and use that
  // as the one material in a LinearMaterial.

  KineticsMaterial::SP_material
    kinmat(new KineticsMaterial(1, 1, 0, "base material"));
  kinmat->set_sigma_t(0, 0,    1.0);
  kinmat->set_sigma_s(0, 0, 0, 1.0);
  LinearMaterial::vec_material materials(1, kinmat);
  LinearMaterial::vec_dbl      times(1, 0.0);
  LinearMaterial::SP_material linmat(new LinearMaterial(times, materials));

  //-------------------------------------------------------------------------//
  // MESH
  //-------------------------------------------------------------------------//

  vec_int fm(1, 10);
  vec_dbl cm(2, 0.0); cm[1] = 10;
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
  for (int i = 0; i < mesh->number_cells()/2; ++i)
    source_map[i] = 1;
  IsotropicSource::spectra_type spectra(2, vec_dbl(1, 0.0));
  // "on"
  spectra[1][0] = 1.0;
  IsotropicSource::SP_externalsource
    q_e1(new IsotropicSource(1, mesh, spectra, source_map));
  // "off"
  spectra[1][0] = 0.0;
  IsotropicSource::SP_externalsource
    q_e2(new IsotropicSource(1, mesh, spectra, source_map));
  //
  // Linear source
  vec_dbl source_times(2, 1.0); source_times[1] = 1.1;
  LinearExternalSource::vec_source sources(2);
  sources.push_back(q_e1);
  sources.push_back(q_e2);
  LinearExternalSource::SP_tdsource
    q_td(new LinearExternalSource(1, mesh, source_times, sources));

  //-------------------------------------------------------------------------//
  // TIME STEPPER
  //-------------------------------------------------------------------------//

  typedef TimeStepper<_1D>            TS_T;


  TS_T stepper(inp, linmat, mesh, false);

  stepper.solve();

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TimeStepper.cc
//---------------------------------------------------------------------------//
