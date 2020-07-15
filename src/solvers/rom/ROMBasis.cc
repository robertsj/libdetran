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
ROMBasis::ROMBasis(int a)
{

}

void ROMBasis::GetBasis(std::string fname, SP_matrix U)
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
    double B[d_n][d_r];
    infile.seekg(0);
    infile.read((char *) &B, sizeof(B)); // read the number of element

    vector<vector<double>> basis_vecs;

    for (int i=0; i<d_r; i++)
    {
      vector<double> v;
      for (int j=0; j< d_n; j++)
	   {
         U->insert(j, i, B[j][i]);
       }
    }
  }
}
}

