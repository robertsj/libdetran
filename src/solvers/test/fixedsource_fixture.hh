//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  fixedsource_fixture.hh
 *  @brief Various data for a fixed source problem
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef FIXEDSOURCE_FIXTURE_HH_
#define FIXEDSOURCE_FIXTURE_HH_

#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "material/Material.hh"
#include "utilities/InputDB.hh"
#include "utilities/MathUtilities.hh"
#include "external_source/ConstantSource.hh"
#include "material/test/material_fixture.hh"

using namespace detran_utilities;
using namespace detran_material;
using namespace detran_geometry;
using namespace detran_external_source;
using namespace detran;

namespace detran_test
{

/// Container for fixed source problem input
struct FixedSourceData
{
  InputDB::SP_input input;
  Material::SP_material material;
  Mesh::SP_mesh mesh;
  ConstantSource::SP_externalsource source;
};

/// Data based on dimension and number of groups
FixedSourceData get_fixedsource_data(unsigned int dim,
                                     unsigned int ng,
                                     unsigned int nfm = 5,
                                     double       cmb = 5.0)
{
  Require(dim >= 0 && dim <= 3);
  Require(ng == 1 || ng == 2 || ng == 7);

  FixedSourceData data;

  // input
  data.input = std::make_shared<InputDB>();
  data.input->put<int>("number_groups",     ng);
  data.input->put<int>("outer_print_level", 1);
  data.input->put<int>("inner_print_level", 0);

  // material
  if (ng == 1)
    data.material = material_fixture_1g();
  else if (ng == 2)
    data.material = material_fixture_2g();
  else if (ng == 7)
    data.material = material_fixture_7g();
  data.material->compute_diff_coef();

  // mesh
  vec_dbl cm(2, 0.0); cm[1] = cmb;
  vec_int fm(1, nfm);
  vec_int mt(1, 0);
  if (dim == 1)
    data.mesh = std::make_shared<Mesh1D>(fm, cm, mt);
  else if (dim == 2)
    data.mesh = std::make_shared<Mesh2D>(fm, fm, cm, cm, mt);
  else if (dim == 3)
    data.mesh = std::make_shared<Mesh3D>(fm, fm, fm, cm, cm, cm, mt);

  // fixed source
  data.source = std::make_shared<ConstantSource>(ng, data.mesh, 1.0);

  return data;
}

/// Return a 7 group slab problem
FixedSourceData get_fixedsource_data_het(unsigned int dim,
                                         unsigned int ng,
                                         unsigned int nfm = 5,
                                         double       cmb = 5.0)
{
  FixedSourceData data;

  // 20 pairs of 1.0 cm slabs of mod and fuel
  vec_dbl cm = detran_utilities::linspace(0, 20, 22);


  data.material = material_fixture_7g();

  return data;
}

} // end namespace detran_test

#endif /* FIXEDSOURCE_FIXTURE_HH_ */

//----------------------------------------------------------------------------//
//              end of fixedsource_fixture.hh
//----------------------------------------------------------------------------//
