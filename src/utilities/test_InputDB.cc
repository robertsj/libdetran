//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_InputDB.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

#include "InputDB.hh"

#include <iostream>
#include <string>

using namespace detran_utils;
using namespace std;

int main()
{
  // Make input
  InputDB::SP_input input;
  input = new InputDB();
  const string s = "number groups";
  input->put(s, 2);
  cout << " number groups = " << input->get<int>(s) << endl;

}

//---------------------------------------------------------------------------//
//              end of test_InputDB.cc
//---------------------------------------------------------------------------//
