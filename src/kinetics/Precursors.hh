//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Precursors.hh
 *  @brief  Precursors
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_PRECURSORS_HH_
#define detran_PRECURSORS_HH_

#include "kinetics/kinetics_export.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

namespace detran
{

/**
 *  @class Precursors
 *  @brief Container for precursor concentrations
 */
class KINETICS_EXPORT Precursors
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Precursors>    SP_precursors;
  typedef detran_utilities::size_t            size_t;
  typedef detran_utilities::vec_dbl           vec_dbl;
  typedef detran_utilities::vec2_dbl          vec2_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_precursor_groups    Number of precursor groups
   *  @param number_cells               Number of cells
   */
  Precursors(const size_t number_precursor_groups,
             const size_t number_cells);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Const accessor to a group precursor field.
   *  @param    i   Group of precursor requested.
   *  @return       Constant reference to precursor
   */
  const vec_dbl& C(const size_t i) const;

  /**
   *  @brief Mutable accessor to a group moments field.
   *  @param    i   Group of precursor requested.
   *  @return       Mutable reference to group moment vector.
   */
  vec_dbl& C(const size_t i);

  size_t number_precursor_groups() const
  {
    return d_number_precursor_groups;
  }

  size_t number_cells() const
  {
    return d_number_cells;
  }

  void display() const;

private:

  /// Number of precursor groups
  size_t d_number_precursor_groups;
  /// Number of spatial cells
  size_t d_number_cells;
  /// Precursor concentrations [npc][ncells]
  vec2_dbl d_C;

};

KINETICS_TEMPLATE_EXPORT(detran_utilities::SP<Precursors>)

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITINS
//---------------------------------------------------------------------------/

#include "Precursors.i.hh"

#endif // detran_PRECURSORS_HH_

//---------------------------------------------------------------------------//
//              end of file Precursors.hh
//---------------------------------------------------------------------------//
