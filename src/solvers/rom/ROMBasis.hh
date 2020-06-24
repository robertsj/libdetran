/*
 * basis.hh
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */

#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <iostream>
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

 // typedef detran_utilities::vec2_dbl     vec2_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

   ROMBasis(const int n, const int r, std::string fname);

  /// virtual destructor
  virtual ~ROMBasis() {};

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  detran_utilities::vec2_dbl GetBasis();

  int rank()
  {
   return d_r;
  }

  int problem_size()
  {
   return d_n;
  }

  //--------------------------------------------------------------------------//
  // PRIVATE MEMBERS
  //--------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
private:
  const int d_r;
  const int d_n;
  std::string d_fname;

};
}
