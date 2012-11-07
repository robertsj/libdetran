//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TransportManager.hh
 *  @brief  TransportManager
 *  @author Jeremy Roberts
 *  @date   Nov 6, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_TRANSPORTMANAGER_HH_
#define detran_TRANSPORTMANAGER_HH_

#include "callow/utils/Initialization.hh"

namespace detran
{

/**
 *  @class TransportManager
 *  @brief Base class for all transport managers
 */
class TransportManager
{

public:

  /// Constructor.  This initializes all external libraries.
  TransportManager(int argc, char** argv)
  {
    callow_initialize(argc, argv);
  }

  /// Destructor.  This initializes all external libraries.
  virtual ~TransportManager()
  {
    callow_finalize();
  }

};


} // end namespace detran

#endif // TRANSPORTMANAGER_HH_ 

//---------------------------------------------------------------------------//
//              end of file TransportManager.hh
//---------------------------------------------------------------------------//
