//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ReactionRates.hh
 * \author robertsj
 * \date   May 24, 2012
 * \brief  ReactionRates class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef REACTIONRATES_HH_
#define REACTIONRATES_HH_

// Detran
#include "FissionSource.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

// System
#include <string>

namespace detran
{

/*!
 *  \class ReactionRates
 *  \brief Computes various reaction rates.
 *
 */
class ReactionRates : public Object
{

public:

  /// \name Useful Typedefs
  /// \{
  typedef SP<ReactionRates>         SP_reactionrates;
  typedef Material::SP_material     SP_material;
  typedef Mesh::SP_mesh             SP_mesh;
  typedef State::SP_state           SP_state;
  typedef FissionSource::SP_source  SP_fissionsource;
  /// \}

  /*!
   *  \brief Constructor.
   *
   *  \param    material    Pin cell pitch (assumed square)
   *  \param    mesh        Vector of fuel pin radii (can be zero length)
   *  \param    state       Region material map (cell-center outward)
   */
  ReactionRates(SP_material material, SP_mesh mesh, SP_state state);

  /// Virtual destructor
  virtual ~ReactionRates(){}

  /// SP Constructor
  static SP<ReactionRates> Create(SP_material material, SP_mesh mesh, SP_state state)
  {
    SP_reactionrates p;
    p = new ReactionRates(material, mesh, state);
    return p;
  }

  /// Set a fission source.
  void set_fission_source(SP_fissionsource fissionsource)
  {
    Require(fissionsource);
    b_fissionsource = fissionsource;
  }

  /*!
   *  \brief Compute pin powers.
   *  \param scale  Normalization for total power (default: unity)
   */
  vec_dbl pin_power(double scale = 1.0);

  /*!
   *  \brief Compute assembly powers.
   *  \param scale  Normalization for total power (default: unity)
   */
  vec_dbl assembly_power(double scale = 1.0);

  /// Verify state correctness.
  bool is_valid() const
  {
    return true;
  }

private:

  /// \name Data
  /// \{

  SP_material b_material;
  SP_mesh     b_mesh;
  SP_state    b_state;
  SP_fissionsource b_fissionsource;

  /// \}


  /// \name Implementation
  /// \{

  /*!
   *  \brief  Generic power function.
   *
   *  \param    key     Mesh map key designating region over which to integrate
   *  \param    scale   Value to which total power is normalized
   */
  vec_dbl power(std::string key, double scale);

  /// \}

};


} // end namespace detran

#endif /* REACTIONRATES_HH_ */
