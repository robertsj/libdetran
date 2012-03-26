//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_State.cc
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
#include "Mesh.hh"
#include "Mesh2D.hh"
#include "Quadrature.hh"
#include "QuadrupleRange.hh"
#include "State.hh"

#include <iostream>
#include <string>

using namespace detran_utils;
using namespace detran;
using namespace std;

int main()
{
  // Make input
  InputDB::SP_input input;
  input = new InputDB();
  input->put<int>("number_groups", 2);
  cout << " number groups = " << input->get<int>("number_groups") << endl;

  // Make geometry.
  Mesh::SP_mesh mesh;
  Mesh::vec_dbl fm(3, 0.0);
  fm[1] = 1.0;
  fm[2] = 2.0;
  Mesh::vec_int mt(4, 1);
  mesh = new Mesh2D(fm, fm, mt);
  cout << " number cells = " << mesh->number_cells() << endl;

  // Make quadrature
  Quadrature::SP_quadrature quadrature;
  quadrature = new QuadrupleRange(2);
  quadrature->display();

  // Make state.
  State::SP_state state;
  state = new State(input, mesh, quadrature);

}

//---------------------------------------------------------------------------//
//              end of test_State.cc
//---------------------------------------------------------------------------//
