//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Homogenize.hh
 *  @brief Homogenize class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_HOMOGENIZE_HH_
#define detran_HOMOGENIZE_HH_

#include "transport/transport_export.hh"
#include "transport/CoarseMesh.hh"
#include "transport/State.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/Definitions.hh"

namespace detran
{

/**
 *  @class Homogenize
 *  @brief Condenses materials on a coarser space and/or energy mesh
 *
 *  Condensation onto coarser space/energy meshes requires a weighting
 *  spectrum.  Two options are provided here.  The first lets the client
 *  use a fine mesh state vector, perhaps from a partially converged
 *  iteration.  The second lets a client define a spectrum for each coarse
 *  mesh region identified by a key.  This allows a user to use e.g. infinite
 *  medium spectra for fine group materials or perhaps homogeneous B0 theory
 *  for assemblies.
 *
 *  Homogenization for the diffusion coefficients is based on flux weighting,
 *  or, if requested and available, current weighting.  Currently, only the
 *  state-based homogenization allows the latter.
 */
/**
 *  @example transport/test/test_Homogenization.cc
 *
 *  Test of Homogenize.
 */
class TRANSPORT_EXPORT Homogenize
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum options_dc
  {
    PHI_D,         // flux-weighted diffusion coefficient
    CURRENT_D,     // current-weighted diffusion coefficient
    PHI_SIGMA_TR,  // flux-weighted transport cross section, D = 1/(3*S_TR)
    END_OPTIONS_DC
  };

  enum options_spectrum
  {
    FINE_MESH_SPECTRUM,   // weight based on fine mesh/group flux
    REGION_SPECTRUM,      // weight based on region mesh/group flux
    END_OPTIONS_SPECTRUM
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef State::SP_state                           SP_state;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef State::vec_moments_type                   vec_moments_type;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::vec2_dbl                vec2_dbl;
  typedef detran_utilities::vec_size_t              vec_size_t;
  typedef detran_utilities::size_t                  size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param material    Original fine group material
   *  @param dc_weight   Diffusion coefficient weighting options
   */
  Homogenize(SP_material material, const size_t dc_weight = 0);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Homogenize the material on a coarser space and energy mesh
   *  @param state          State vector for flux weighting
   *  @param mesh           Fine mesh with appropriate coarse mesh map
   *  @param regionkey      Edit region mesh map key
   *  @param coarsegroup    Optional vector of fine groups per coarse group
   *  @param groups         Optional vector of groups over which to homogenize
   */
  SP_material homogenize(SP_state           state,
                         SP_mesh            mesh,
                         const std::string &regionkey,
                         vec_int            coarsegroup = vec_int(0),
                         vec_size_t         groups = vec_size_t(0));

  /**
   *  @brief Homogenize the material on a coarser space and energy mesh
   *  @param spectrum       Region-wise spectrum for flux weighting [nr][ng]
   *  @param spectrumkey    Spectrum region mesh map key
   *  @param mesh           Fine mesh with appropriate mesh maps
   *  @param regionkey      Edit region mesh map key
   *  @param coarsegroup    Optional vector of fine groups per coarse group
   *  @param groups         Optional vector of groups over which to homogenize
   */
  SP_material homogenize(const vec2_dbl    &spectrum,
                         const std::string &spectrumkey,
                         SP_mesh            mesh,
                         const std::string &regionkey,
                         vec_int            coarsegroup = vec_int(0),
                         vec_size_t         groups = vec_size_t(0));

  /// Reset the weighting method for diffusion coefficient generation
  void set_option_dc(const size_t option_dc);

  /// Get the coarse mesh flux
  const vec2_dbl& coarse_mesh_flux(){return d_coarse_mesh_flux;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Original fine group material
  SP_material d_material;
  /// Original number of groups
  size_t d_number_groups;
  /// Flag to specify method for diffusion coefficient generation
  size_t d_option_dc;
  /// Flag to specify method for mesh/group condensation
  size_t d_option_spectrum;
  /// State vector for weighting
  SP_state d_state;
  /// Region spectrum for weighting [ng][nr]
  vec2_dbl d_spectrum;
  /// Fine mesh spectrum map
  vec_int d_spectrum_map;
  /// Coarse mesh flux
  vec2_dbl d_coarse_mesh_flux;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /**
   *  @brief Homogenize a material
   *  @param mesh           Fine mesh with appropriate coarse mesh map
   *  @param key            Coarse mesh map key
   *  @param coarsegroup    Vector of fine groups per coarse group
   */
  SP_material homogenize(SP_mesh            mesh,
                         const std::string &region,
                         vec_int            coarsegroup,
                         vec_size_t         groups);

  /**
   *  @brief Return the flux spectrum used for weighting in a cell and group
   *  @param cell     Mesh cell being evaluated
   *  @param g        Energy group being evaluated
   */
  double spectrum(const size_t g, const size_t cell) const;

  /**
   *  @brief Return the current spectrum used for weighting in a cell and group
   *  @param cell     Mesh cell being evaluated
   *  @param g        Energy group being evaluated
   */
  double current(const size_t g, const size_t cell) const;

};

} // end namespace detran

#endif // detran_HOMOGENIZE_HH_

//----------------------------------------------------------------------------//
//              end of file Homogenize.hh
//----------------------------------------------------------------------------//
