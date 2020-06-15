/*
 * basis.hh
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */

#ifndef SOLVERS_ROM_BASIS_HH_
#define SOLVERS_ROM_BASIS_HH_


#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "utilities/Definitions.hh"

#include <iostream>
using namespace std;


namespace detran
{
//----------------------------------------------------------------------------//
/**
 *  @class basis
 *  @brief get the basis from a binary file and store them in a matrix
 */
//----------------------------------------------------------------------------//
class Basis
{

public:
  typedef detran_utilities::vec2_dbl     vec2_dbl;

  // constructor
  Basis(const int d_n, const int d_r, string fname);

  vec2_dbl GetBasis(std::string fname, int r);

protected:

  //rank
  const int d_r;
  const int d_n;
  std::string fname;

};






















}











#endif /* SOLVERS_ROM_BASIS_HH_ */
