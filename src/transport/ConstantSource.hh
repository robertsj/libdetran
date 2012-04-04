//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConstantSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  ConstantSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef CONSTANTSOURCE_HH_
#define CONSTANTSOURCE_HH_

#include "ExternalSource.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class ConstantSource
 * \brief A single isotropic source everywhere.
 */
//===========================================================================//

class ConstantSource: public ExternalSource
{

public:

  // Source types
  typedef detran_utils::SP<ConstantSource>  SP_source;

  IsotropicSource(SP_mesh mesh,
                  SP_quadrature quadrature,
                  int number_groups)
    :  ExternalSource(mesh, quadrature, number_groups)
  { /* ... */ }

  virtual double source(int cell, int group)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    return d_source;
  }

  virtual double source(int cell, int group, int angle)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    Require(angle >=0);
    Require(angle < d_number_angles);
    return d_source * d_quadrature->angular_norm();
  }

  void set_source(double source)
  {
    d_source = source;
  }

private:

  double d_source;

};

} // end namespace detran

#endif /* CONSTANTSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of ConstantSource.hh
//---------------------------------------------------------------------------//
