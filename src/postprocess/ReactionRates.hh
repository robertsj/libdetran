//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ReactionRates.hh
 *  @author robertsj
 *  @date   May 24, 2012
 *  @brief  ReactionRates class definition.
 */
//---------------------------------------------------------------------------//

#ifndef REACTIONRATES_HH_
#define REACTIONRATES_HH_

#include "postprocess/postprocess_export.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/FissionSource.hh"
#include "transport/State.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <string>

namespace detran_postprocess
{

/**
 *  @class ReactionRates
 *  @brief Computes various reaction rates based on the state.
 *
 *  Once the problem is solved, several quantities are often
 *  required for analysis.  These include global net gains and
 *  losses, reaction rates in edit regions such as pins or
 *  assemblies, and so on.
 *
 */
class POSTPROCESS_EXPORT ReactionRates
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<ReactionRates>       SP_reactionrates;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran::State::SP_state                   SP_state;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::vec_int                 vec_int;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *
   *  @param    material    Pin cell pitch (assumed square)
   *  @param    mesh        Vector of fuel pin radii (can be zero length)
   *  @param    state       Region material map (cell-center outward)
   */
  ReactionRates(SP_material material, SP_mesh mesh, SP_state state);

  /// Virtual destructor
  virtual ~ReactionRates(){}

  /// SP Constructor
  static SP_reactionrates
  Create(SP_material material, SP_mesh mesh, SP_state state)
  {
    SP_reactionrates p;
    p = new ReactionRates(material, mesh, state);
    return p;
  }

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Relative power of an edit region.
   *
   *  Note, the region powers are returned as a vector
   *  in the same order as the regions are indexed.  That means
   *  the user is responsible for mapping the region powers
   *  to the end application.
   *
   *  For the built-in pin and assembly arrays, the indexing
   *  is natural, following the same x then y then z ordering
   *  as used in \ref Mesh.
   *
   *  @param key    String identifier for the edit region
   *  @param scale  Total power for normalization.  Negative for no scaling.
   */
  vec_dbl region_power(std::string key, double scale = 1.0);

  /**
   *  @brief Edit mesh function
   *
   *  Given a function defined on the fine mesh, return the
   *  function defined within an edit regions defined by the key.
   *  The user can optionally average the value (rather than
   *  just integrate it)
   *
   *  @param key                String identifier for the edit region
   *  @param fine_mesh_function Fine mesh function
   *  @param mean               Compute the mean in the edit region
   */
  vec_dbl edit(std::string key,
               const vec_dbl &fine_mesh_function,
               bool mean = false);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  SP_material b_material;
  SP_mesh     b_mesh;
  SP_state    b_state;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

};


} // end namespace detran_postprocess

#endif /* REACTIONRATES_HH_ */
