//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Homogenize.hh
 *  @brief  Homogenize 
 *  @author Jeremy Roberts
 *  @date   Mar 26, 2013
 */
//---------------------------------------------------------------------------//

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
 *  Homogenization is based on flux or current weighting.
 */
class TRANSPORT_EXPORT Homogenize
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum diff_coef_weighting
  {
    PHI_D,         // flux-weighted diffusion coefficient
    CURRENT_D,     // current-weighted diffusion coefficient
    PHI_SIGMA_TR,  // flux-weighted transport cross section, D = 1/(3*S_TR)
    END_DIFF_COEF_WEIGHTING
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef State::SP_state                           SP_state;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::vec2_dbl                vec2_dbl;
  typedef detran_utilities::size_t                  size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param material    Original fine group material
   */
  Homogenize(SP_material material);

  /**
   *  @brief Homogenize the material on a coarser space and energy mesh
   *  @param state          State vector for flux weighting
   *  @param mesh           Fine mesh with appropriate coarse mesh map
   *  @param key            Coarse mesh map key
   *  @param coarsegroup    Vector of fine groups per coarse group
   *  @param dc_weight      Selects weighting scheme for diffusion coefficient
   */
  SP_material homogenize(SP_state state,
                         SP_mesh mesh,
                         std::string key,
                         vec_int coarsegroup,
                         size_t dc_weight = 0);

  /**
   *  @brief Homogenize the material on a coarser space mesh.
   *  @param state          State vector for flux weighting
   *  @param mesh           Fine mesh with appropriate coarse mesh map
   *  @param key            Coarse mesh map key
   *  @param dc_weight      Selects weighting scheme for diffusion coefficient
   */
  SP_material homogenize(SP_state state,
                         SP_mesh mesh,
                         std::string key,
                         size_t dc_weight = 0);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Original fine group material
  SP_material d_material;
  /// Original number of groups
  size_t d_number_groups;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Return the current (or flux) vector from state
  const vec_dbl& current(SP_state, size_t g, size_t dc_weight = 0) const;

};

} // end namespace detran

#endif // detran_HOMOGENIZE_HH_

//---------------------------------------------------------------------------//
//              end of file Homogenize.hh
//---------------------------------------------------------------------------//
