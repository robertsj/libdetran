/*
 * basis.cc
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */

#include "ROMBasis.hh"
#include <fstream>

using namespace std;

namespace detran
{
ROMBasis::ROMBasis()
{

}

void ROMBasis::GetBasis(const char* fname, SP_matrix U)
{
 int d_r = U->number_columns();
 int d_n = U->number_rows();

 ifstream infile;
 infile.open(fname, ios::binary | ios::in);

 if(!infile)
 {
  cout << "Cannot open file!" << endl;
 }
 else
 {
  double * B = new double[d_r*d_n];
  infile.seekg(0);
  infile.read((char *) B, (d_n*d_r)*sizeof(double)); // read the number of element
  int i;
  int j;
  for (int c=0; c<d_r*d_n; c++)
  {
   i = c/d_r;
   j = c%d_r;
   U->insert(i, j, B[c]);
   }
 }
}
}

