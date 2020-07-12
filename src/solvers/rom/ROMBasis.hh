/*
 * basis.hh
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */

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

  ROMBasis(int a);

  /// virtual destructor
  virtual ~ROMBasis() {};

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  void GetBasis(std::string fname, SP_matrix U);

  /*
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
  */

};
}

#endif //  detran_rombasis_HH_
