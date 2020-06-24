/*
 * basis.cc
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */



#include <ROMBasis.hh>
#include <fstream>

using namespace std;

namespace detran
{
ROMBasis::ROMBasis(const int n, const int r, string fname)
:d_n(n),
 d_r(r),
 d_fname(fname)
{
  std::cout << "******** please work ********" << "\n";
}


detran_utilities::vec2_dbl ROMBasis::GetBasis()
{
  ifstream infile;
  infile.open(d_fname, ios::binary | ios::in);

  double U[d_n][d_r];
  infile.seekg(0);
  infile.read((char *) &U, sizeof(U)); // read the number of element

  vector<vector<double>> basis_vecs;

  for (int i=0; i<d_r; i++)
  {
    vector<double> v;
    for (int j=0; j< d_n; j++)
	{
	 //std::cout << U[j][i] << "\n";
	 v.push_back(U[j][i]);
    }
	 basis_vecs.push_back(v);
  }

  return  basis_vecs;
}

}


