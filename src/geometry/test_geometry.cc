//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_geometry.cc
 * \author Jeremy Roberts
 * \date   Mar 20, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#include <vector>
#include <iostream>
#include "Mesh2D.hh"

using namespace std;
using namespace detran;

int main()
{

  vector<int>     xfm(10, 2);
  vector<int>     yfm(10, 2);
  vector<double>  xcm(11, 0.0);
  vector<double>  ycm(11, 0.0);
  vector<int>     mat_map(100, 0);
  Mesh2D mesh(xfm, yfm, xcm, ycm, mat_map);

  cout << " dx(1) = " << mesh.dx(1) << endl;

}


//---------------------------------------------------------------------------//
//              end of test_geometry.cc
//---------------------------------------------------------------------------//
