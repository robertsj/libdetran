//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ExternalSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  ExternalSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EXTERNALSOURCE_HH_
#define EXTERNALSOURCE_HH_

// Detran
#include "Mesh.hh"
#include "Quadrature.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class ExternalSource
 * \brief Base volume source class.
 */
//===========================================================================//

class ExternalSource: public Object
{

public:

  typedef SP<ExternalSource>        SP_source;
  typedef Mesh::SP_mesh             SP_mesh;
  typedef Quadrature::SP_quadrature SP_quadrature;


  ExternalSource(SP_mesh mesh,
                 SP_quadrature quadrature,
                 int number_groups)
    : d_mesh(mesh)
    , d_quadrature(quadrature)
    , d_number_groups(number_groups)
    , d_number_angles(quadrature->number_angles())
  {}

  /// Virtual destructor
  virtual ~ExternalSource(){}

  /*!
   *  \brief Get moments source for cell.
   *
   */
  virtual double source(int cell, int group) = 0;

  /*!
   *  \brief Get discrete source for cell and cardinal angle.
   *
   */
  virtual double source(int cell, int group, int angle) = 0;

  virtual bool is_valid() const
  {
    /* ... */
  }

protected:

  /// Cartesian mesh.
  SP_mesh d_mesh;

  /// Quadrature
  SP_quadrature d_quadrature;

  /// Number of groups.
  const int d_number_groups;

  /// Number of angles
  const int d_number_angles;

  /// Am I ready?
  int d_initialized;

};

} // end namespace detran

#endif /* EXTERNALSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of ExternalSource.hh
//---------------------------------------------------------------------------//
