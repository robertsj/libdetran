//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumPinCell.cc
 *  @brief SpectrumPinCell member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "SpectrumPinCell.hh"
#include "EigenvalueManager.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/solver/Eispack.hh"
#include "utilities/MathUtilities.hh"

namespace detran
{

//----------------------------------------------------------------------------//
SpectrumPinCell::SpectrumPinCell(SP_input    input,
                                 SP_material material,
                                 SP_mesh     mesh)
  : SpectrumBase(input, material, mesh)
{
  Insist(d_input->check("mgpc_spectrum_pincell_db"),
         "Pincell spectrum requires user-defined pincells for each region");

  std::string key = "PINCELL_SPECTRUM";

  // Get the pincell db and set things that must be set from other data
  SP_input db = d_input->get<SP_input>("mgpc_spectrum_pincell_db");
  db->put<int>("number_groups", d_material->number_groups());
  db->put<std::string>("bc_west",  "reflect");
  db->put<std::string>("bc_east",  "reflect");
  db->put<std::string>("bc_south", "reflect");
  db->put<std::string>("bc_north", "reflect");
  db->put<int>("inner_print_level", 0);
  db->put<int>("outer_print_level", 0);


  // Determine the number of regions with this key in the mesh map
  Insist(mesh->mesh_map_exists(key), "Mesh map must exist for spectrum");
  const vec_int &mesh_map = mesh->mesh_map(key);
  vec_int unique = detran_utilities::vec_unique(mesh_map);
  size_t npins = unique.size();

  // Define the pincells
  d_pincells.resize(npins);
  for (size_t p = 0; p < npins; ++p)
  {
    Assert(db->check("pitch_"+AsString(p)));
    double w = db->get<double>("pitch_"+AsString(p));
    detran_geometry::Point pitch(w, w);

    Assert(db->check("mat_map_"+AsString(p)));
    vec_int mat  = db->get<vec_int>("mat_map_"+AsString(p));

    Assert(db->check("radii_"+AsString(p)));
    vec_dbl rad  = db->get<vec_dbl>("radii_"+AsString(p));

    d_pincells[p] = detran_geometry::PinCell::Create(pitch, mat, rad);
    d_pincells[p]->meshify(20, false);
  }

  d_pincell_db = db;
}

//----------------------------------------------------------------------------//
SpectrumPinCell::vec2_dbl SpectrumPinCell::spectrum(double keff)
{
  using namespace callow;
  using namespace detran_material;

  size_t ng = d_material->number_groups();
  size_t np = d_pincells.size();

  vec2_dbl xi(ng, vec_dbl(np, 0.0));

  Material::SP_material mat(new Material(*d_material));
  for (size_t m = 0; m < mat->number_materials(); ++m)
  {
    for (size_t g = 0; g < ng; ++g)
    {
      double sf = d_material->nu_sigma_f(m, g) / keff;
      mat->set_sigma_f(m, g, sf);
    }
  }

  for (size_t p = 0; p < np; ++p)
  {
    SP_mesh pinmesh = d_pincells[p]->mesh();
    EigenvalueManager<_2D> M(d_pincell_db, mat, pinmesh);
    M.solve();
    callow::Vector xi_p(ng, 0.0);
    for (size_t g = 0; g < ng; ++g)
    {
      for (size_t i = 0; i < pinmesh->number_cells(); ++i)
      {
        xi_p[g] += M.state()->phi(g)[i] * pinmesh->volume(i);
      }
    }
    xi_p.scale(1.0 / xi_p.norm());
    xi_p.display();
    for (size_t g = 0; g < ng; ++g)
      xi[g][p] = xi_p[g];
  }
  return xi;
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file SpectrumPinCell.cc
//----------------------------------------------------------------------------//

