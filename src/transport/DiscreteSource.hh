//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiscreteSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  DiscreteSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef DISCRETESOURCE_HH_
#define DISCRETESOURCE_HH_

#include "ExternalSource.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class DiscreteSource
 * \brief 
 */
//===========================================================================//

class DiscreteSource: public ExternalSource
{

public:

  // Source types
  typedef SP<DiscreteSource>  SP_source;

  DiscreteSource(SP_mesh        mesh,
                 SP_quadrature  quadrature,
                 int            number_groups,
                 vec3_dbl       q)
    : ExternalSource(mesh, quadrature, number_groups)
    , d_source(q)
  {
    Require(d_source.size()       == number_groups);
    Require(d_source[0].size()    == quadrature->number_angles());
    Require(d_source[0][0].size() == mesh->number_cells());

  }

  static SP<ExternalSource>
  Create(SP<detran::Mesh>       mesh,
         SP<detran::Quadrature> quadrature,
         int                    number_groups,
         vec3_dbl               q)
  {
    SP_source p;
    p = new DiscreteSource(mesh, quadrature, number_groups, q);
    return p;
  }

  virtual double source(int cell, int group)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    double value = 0;
    for (int o = 0; o < d_quadrature->number_octants(); o++)
    {
      for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
      {
        int angle = d_quadrature->index(o, a);
        value += d_source[group][angle][cell] * d_quadrature->weight(a);
      }
    }
    return value;
  }

  virtual double source(int cell, int group, int angle)
  {
    Require(cell >= 0);
    Require(cell < d_mesh->number_cells());
    Require(group >= 0);
    Require(group < d_number_groups);
    Require(angle >=0);
    Require(angle < d_number_angles);
    return d_source[group][angle][cell];
  }

//  void set_source(vec3_dbl &source)
//  {
//    Require(source.size() == d_number_groups);
//    Require(source[0].size() == d_number_angles);
//    Require(source[0][0].size() == d_mesh->number_cells());
//    d_source = source;
//  }

private:

  /// Source for all points, angles, and energies
  vec3_dbl d_source;

};

} // end namespace detran

#endif /* DISCRETESOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of DiscreteSource.hh
//---------------------------------------------------------------------------//
