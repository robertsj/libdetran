//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  eigenvalue_fixture.hh
 *  @brief Various data for an eigenvalue problem
 */
//----------------------------------------------------------------------------//

#ifndef EIGENVALUE_FIXTURE_HH_
#define EIGENVALUE_FIXTURE_HH_

#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "material/Material.hh"
#include "utilities/InputDB.hh"
#include "material/test/material_fixture.hh"

using namespace detran_utilities;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran;

namespace detran_test
{

/// Container for fixed source problem input
struct EigenvalueData
{
  InputDB::SP_input input;
  Material::SP_material material;
  Mesh::SP_mesh mesh;
};

/// Data based on dimension and number of groups
EigenvalueData get_eigenvalue_data(unsigned int dim, unsigned int ng)
{
  Require(dim >= 0 && dim <= 3);
  Require(ng == 1 || ng == 2 || ng == 7);

  EigenvalueData data;

  // input
  data.input = InputDB::Create();
  data.input->put<int>("number_groups",      ng);
  data.input->put<int>("outer_print_level",  1);
  data.input->put<int>("inner_print_level",  0);
  data.input->put<std::string>("bc_west",    "reflect");
  data.input->put<std::string>("bc_east",    "reflect");
  data.input->put<double>("inner_tolerance", 1e-18);
  data.input->put<double>("outer_tolerance", 1e-16);
  data.input->put<double>("eigen_tolerance", 1e-16);
  data.input->put<int>("inner_max_iters",    1000000);
  data.input->put<int>("outer_max_iters",    1000000);
  data.input->put<int>("eigen_max_iters",    1000000);

  InputDB::SP_input callow_db = InputDB::Create();
  callow_db->put<std::string>("eigen_solver_type", "gd");
  callow_db->put<int>("eigen_solver_monitor_level", 2);
  callow_db->put<double>("eigen_solver_tol", 1.0e-14);
  data.input->put<InputDB::SP_input>("eigen_solver_db", callow_db);

  // material
  if (ng == 1)
    data.material = material_fixture_1g(); //
  else if (ng == 2)
    data.material = material_fixture_2g();
  else if (ng == 7)
    data.material = material_fixture_7g();
  data.material->compute_diff_coef();

  // mesh
  vec_dbl cm(2, 0.0); cm[1] = 5.0;
  vec_int fm(1, 5);
  vec_int mt(1, 2);
  if (dim == 1)
    data.mesh = Mesh1D::Create(fm, cm, mt);
  else if (dim == 2)
    data.mesh = Mesh2D::Create(fm, fm, cm, cm, mt);
  else if (dim == 3)
    data.mesh = Mesh3D::Create(fm, fm, fm, cm, cm, cm, mt);

  return data;
}

} // end namespace detran_test

#endif /* EIGENVALUE_FIXTURE_HH_ */

//----------------------------------------------------------------------------//
//              end of eigenvalue_fixture.hh
//----------------------------------------------------------------------------//
