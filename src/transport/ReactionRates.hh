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
#include "Boundary.hh"
#include "FissionSource.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

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
    p = new ReactionRates(SP_material material, SP_mesh mesh, SP_state state);
    return p;
  }

  /// Verify state correctness.
  bool is_valid()
  {
    return true;
  }

};


} // end namespace detran

#endif /* REACTIONRATES_HH_ */
