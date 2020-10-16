//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  OrthogonalBasisParameters.hh
 *  @brief OrthogonalBasisParameters struct definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_orthog_ORTHOGONALBASISPARAMETERS_HH_
#define detran_orthog_ORTHOGONALBASISPARAMETERS_HH_

#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"
#include <string>

namespace detran_orthog
{

// Convenience structure for passing parameters to the basis
class OrthogonalBasisParameters
{
public:
  typedef detran_utilities::size_t                        size_t;
  typedef detran_utilities::vec_dbl                       vec_dbl;

  OrthogonalBasisParameters()
    : order(0)
    , size(0)
    , orthonormal(false)
    , even_only(false)
    , lower_bound(-1)
    , upper_bound(1)
    , transformed_key("dct")
    , transformed_option(0)
    , db(new detran_utilities::InputDB())
  {}
  /// Order of basis
  size_t order;
  /// Size of basis vector
  size_t size;
  /// Orthonormalize
  bool orthonormal;
  /// Use even orders only
  bool even_only;
  /// Points and weights for continuous bases
  //@{
  detran_utilities::vec_dbl x;
  detran_utilities::vec_dbl qw;
  //@}
  /// Weights for discrete basis
  //@{
  detran_utilities::vec_dbl w;
  //@}
  /// Lower and upper bounds for continuous bases
  //@{
  double lower_bound;
  double upper_bound;
  //@}
  /// Key for transformed basis
  std::string transformed_key;
  /// Option for transformed basis
  size_t transformed_option;

  detran_utilities::InputDB::SP_input db;
};

} // end namespace detran_orthog

#endif /* detran_orthog_ORTHOGONALBASISPARAMETERS_HH_ */

//----------------------------------------------------------------------------//
//              end of file OrthogonalBasisParameters.hh
//----------------------------------------------------------------------------//
