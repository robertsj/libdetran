//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumBase.hh
 *  @brief SpectrumBase class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_SPECTRUMBASE_HH_
#define detran_SPECTRUMBASE_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"

namespace detran
{

/**
 *  @class SpectrumBase
 *  @brief Base class for spectra used in coarse mesh preconditioning
 */
class SpectrumBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::vec2_dbl                vec2_dbl;
  typedef detran_utilities::vec_size_t              groups_t;
  typedef groups_t::iterator                        groups_iter;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  SpectrumBase(SP_input input, SP_material material, SP_mesh mesh)
    : d_input(input)
    , d_material(material)
    , d_mesh(mesh)
    , d_include_fission(false)
  {
    Require(input);
    Require(material);
    Require(mesh);
    if (d_input->check("mgpc_spectrum_include_fission"))
    {
      d_include_fission =
        0 != d_input->get<int>("mgpc_spectrum_include_fission");
    }
  }

  /// Virtual destructor
  virtual ~SpectrumBase(){}

  /// Produce the edit region-dependent spectra
  virtual vec2_dbl spectrum(double keff = 1.0) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_input d_input;
  SP_material d_material;
  SP_mesh d_mesh;
  bool d_include_fission;

};

} // end namespace detran

#endif /* detran_SPECTRUMBASE_HH_ */

//----------------------------------------------------------------------------//
//              end of file SpectrumBase.hh
//----------------------------------------------------------------------------//
