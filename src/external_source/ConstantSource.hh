//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ConstantSource.hh
 *  @brief ConstantSource class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_external_source_CONSTANTSOURCE_HH_
#define detran_external_source_CONSTANTSOURCE_HH_

#include "ExternalSource.hh"

namespace detran_external_source
{

//----------------------------------------------------------------------------//
/**
 *  @class ConstantSource
 *  @brief Defines a single isotropic source everywhere.
 */
//----------------------------------------------------------------------------//

class EXTERNAL_SOURCE_EXPORT ConstantSource: public ExternalSource
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param strength       Source strength in all groups and all space;
   *                        the unit is n/cc-sec
   *  @param quadrature     Pointer to quadrature (optional)
   */
  ConstantSource(size_t number_groups,
                 SP_mesh mesh,
                 double strength,
                 SP_quadrature quadrature = SP_quadrature(0));

  /// SP constructor
  static SP_externalsource
  Create(size_t number_groups,
         SP_mesh mesh,
         double source,
         SP_quadrature quadrature = SP_quadrature(0));

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EXTERNAL SOURCES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  double source(const size_t cell, const size_t group)
  {
    Require(cell < d_mesh->number_cells());
    Require(group < d_number_groups);
    return d_source;
  }

  double source(const size_t cell, const size_t group, const size_t angle)
  {
    Require(cell < d_mesh->number_cells());
    Require(group < d_number_groups);
    Require(angle < d_number_angles);
    return d_discrete_source;
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Source strength (in all groups, in all space)
  const double d_source;

  /// Source strength * angular norm
  const double d_discrete_source;

};

} // end namespace detran_external_source

#endif /* detran_external_source_CONSTANTSOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of ConstantSource.hh
//----------------------------------------------------------------------------//
