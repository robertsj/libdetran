/*
 * basis.cc
 *
 *  Created on: Jun 15, 2020
 *      Author: rabab
 */



#include "utilities/Warning.hh"
#include "Basis.hh"
#include <fstream>

using namespace std;
typedef detran_utilities::vec2_dbl     vec2_dbl;


namespace detran
{

Basis::Basis(const int d_n, const int d_r, string fname)
:d_n(1),
 d_r(1)

{
double U[d_n][d_r];
};

vec2_dbl Basis::GetBasis(std::string fname, int r)
{
  ifstream infile;
  infile.open(fname, ios::binary | ios::in);

  double U[d_n][d_r];
  infile.seekg(0);
  infile.read((char *) &U, sizeof(U)); // read the number of element

  vector<vector<double>> basis_vecs;

  for (int i=0; i<r; i++)
  {
    vector<double> v;
    for (int j=0; j< d_n; j++)
	{
	 std::cout << U[j][i] << "\n";
	 v.push_back(U[j][i]);
     }
	 basis_vecs.push_back(v);
  }

  return  basis_vecs;
}


/*
MatrixDense Basis::ToMatdense()
{

 callow::MatrixDense::SP_matrix A_r;
 A_r = new callow::MatrixDense(5, 5);
 for (int i=0; i<5; i++)
 {
   for (int k=0; k<5 ; k++)
   {
    for (int j=0; j<20; j++)
    {
      double v = U[j][i]*A_r(j, k);
      //std::cout << U[j][i] << " & ";
      A_r->insert(i, k, v, 1);
    }
   }
 }

}

*/


}
