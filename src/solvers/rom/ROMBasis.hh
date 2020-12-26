//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ROMBasis.hh
 *  @brief RoMBasis class definition.
 *  @note  Copyright(C) 2020 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_rombasis_HH_
#define detran_rombasis_HH_


#include "callow/matrix/MatrixDense.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <iostream>
#include <string>
#include <vector>


namespace detran
{
//----------------------------------------------------------------------------//
/**
 *  @class basis
 *  @brief get the basis from a binary file and store them in a matrix
 */
//----------------------------------------------------------------------------//

class ROMBasis
{

public:
  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//
  typedef callow::MatrixDense::SP_matrix             SP_matrix;
  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  ROMBasis();

  // virtual destructor
  virtual ~ROMBasis() {};

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Read the basis from file and store it in a dense matrix
  static void GetBasis(const char* fname, SP_matrix U);

};
}

#endif //  detran_rombasis_HH_
