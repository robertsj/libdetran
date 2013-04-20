//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundarySource.hh
 * \brief  BoundarySource class definition
 * \author Jeremy Roberts
 * \date   Sep 10, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYSOURCE_HH_
#define BOUNDARYSOURCE_HH_

#include "utilities/Definitions.hh"

namespace detran_external_source
{

/*!
 *  \class BoundarySource
 *  \brief Base boundary source class
 */
class BoundarySource
{

public:

  typedef detran_utilities::size_t size_t;

protected:

  /// Index of surface on which the source is defined
  size_t d_surface;

};


} // end namespace detran_external_source

#endif // BOUNDARYSOURCE_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundarySource.hh
//---------------------------------------------------------------------------//
