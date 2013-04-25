//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Execute.hh
 *  @brief  Execute class definition.
 *  @note   Copyright (C) 2013 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/** @mainpage detran: A DEterministic TRANsport package
 *
 *  @section sec_introduction Introduction
 *
 *  This is the introduction.
 *
 *  @section install_sec Installation
 *
 *  @subsection step1 Step 1: Opening the box
 *
 *  etc...
 */
//---------------------------------------------------------------------------//

#ifndef detran_EXECUTE_HH_
#define detran_EXECUTE_HH_

#include "detran_config.hh"
#include "StupidParser.hh"
#include "external_source/ExternalSource.hh"
#include "external_source/ConstantSource.hh"
#include "external_source/DiscreteSource.hh"
#include "external_source/IsotropicSource.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "solvers/FixedSourceManager.hh"
#include "solvers/EigenvalueManager.hh"
#include "transport/DimensionTraits.hh"
#include "transport/State.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include <iostream>
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Execute
 *  @brief Setup and execute the problem.
 */
//---------------------------------------------------------------------------//

class Execute
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input             SP_input;
  typedef detran_geometry::Mesh::SP_mesh                  SP_mesh;
  typedef detran_material::Material::SP_material          SP_material;
  typedef detran_external_source::
          ExternalSource::SP_externalsource               SP_externalsource;
  typedef detran_external_source::
          ExternalSource::vec_externalsource              vec_externalsource;
  typedef detran::State::SP_state                         SP_state;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*
   *  \brief Constructor
   *  \param parser     Input parser
   */
  Execute(StupidParser &parser);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Solve the problem.
  template <class D>
  void solve();

  /// Write output to file.
  void output();

  // Return dimension of problem.
  int dimension()
  {
    return d_dimension;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Input
  SP_input d_input;
  // Material
  SP_material d_material;
  // Mesh
  SP_mesh d_mesh;
  // State vector
  SP_state d_state;
  // Number of groups
  int d_number_groups;
  // External source
  SP_externalsource d_externalsource;
  // Problem dimensions
  int d_dimension;
  // Problem type
  std::string d_problem_type;

};

} // end namespace detran

#endif /* detran_EXECUTE_HH_ */

//---------------------------------------------------------------------------//
//              end of Execute.hh
//---------------------------------------------------------------------------//
