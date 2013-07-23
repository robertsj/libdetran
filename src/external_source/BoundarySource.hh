//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  BoundarySource.hh
 *  @brief BoundarySource class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_external_source_BOUNDARYSOURCE_HH_
#define detran_external_source_BOUNDARYSOURCE_HH_

#include "utilities/Definitions.hh"

namespace detran_external_source
{

/**
 *  @class BoundarySource
 *  @brief Base boundary source class
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

#endif // detran_external_source_BOUNDARYSOURCE_HH_

//----------------------------------------------------------------------------//
//              end of file BoundarySource.hh
//----------------------------------------------------------------------------//
