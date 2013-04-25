//----------------------------------*-C++-*----------------------------------//
/**
 * @file   MultiPhysics.hh
 * @brief  MultiPhysics
 * @author Jeremy Roberts
 * @date   Dec 19, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_MULTIPHYSICS_HH_
#define detran_MULTIPHYSICS_HH_

#include "kinetics/kinetics_export.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <map>

namespace detran
{

/**
 *  @class MultiPhysics
 *  @brief Container for multiphysics state variables
 *
 *  This class wraps a map of double vectors that represent various
 *  state variables.  The class is "dumb" in that it knows nothing
 *  about the problem as is merely a data container.
 *
 */
class KINETICS_EXPORT MultiPhysics
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<MultiPhysics>    SP_multiphysics;
  typedef detran_utilities::size_t              size_t;
  typedef detran_utilities::vec_dbl             vec_dbl;
  typedef detran_utilities::vec2_dbl            vec2_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  MultiPhysics(size_t number_variables = 0);

  /// SP constructor
  static SP_multiphysics Create(size_t number_variables = 0);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Const accessor to a physics variable
   *  @param    id    Identifier
   *  @return         Constant reference to physics variable
   */
  const vec_dbl& variable(const size_t id) const;

  /**
   *  @brief Mutable accessor to a physics variable
   *  @param    id    Identifier
   *  @return         Mutable reference to physics variable
   */
  vec_dbl& variable(const size_t id);

  /**
   *  @brief Add a physics variable
   *  @param    id   String identifier
   *  @param    value Physics value
   */
  void add_variable(const size_t id, vec_dbl value);

  /// Return the number of physics variables we have.
  size_t number_variables() const
  {
    return d_number_variables;
  }

  /// Pretty display of contents
  void display() const;

private:

  /// Number of variables
  size_t d_number_variables;
  /// Variables
  vec2_dbl d_physics_variables;

};

KINETICS_TEMPLATE_EXPORT(detran_utilities::SP<MultiPhysics>)

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITINS
//---------------------------------------------------------------------------/

//#include "MultiPhysics.i.hh"


#endif // detran_MULTIPHYSICS_HH_

//---------------------------------------------------------------------------//
//              end of file MultiPhysics.hh
//---------------------------------------------------------------------------//
