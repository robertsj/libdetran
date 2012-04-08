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
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Quadrature::SP_quadrature         SP_quadrature;

  ConstantSource(SP_mesh mesh,
                 SP_quadrature quadrature,
                 int number_groups)
    :  ExternalSource(mesh, quadrature, number_groups)
  { /* ... */ }

  /*!
   *  \brief SP Constructor.
   */
  static detran_utils::SP<ExternalSource>
  Create(detran_utils::SP<detran::Mesh> mesh,
         detran_utils::SP<detran::Quadrature> quadrature,
         int number_groups,
         double source)
  {
    ConstantSource::SP_source p;
    p = new ConstantSource(mesh, quadrature, number_groups);
    p->set_source(source);
    return p;
  }

  double source(int cell, int group)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    return d_source;
  }

  double source(int cell, int group, int angle)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    Require(angle >=0);
    Require(angle < d_number_angles);
    return d_source * 1.0; //d_quadrature->angular_norm();
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
