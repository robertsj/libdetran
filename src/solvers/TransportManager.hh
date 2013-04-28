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

#include "solvers/solvers_export.hh"
#include "callow/utils/Initialization.hh"

namespace detran
{

/**
 *  @class TransportManager
 *  @brief Base class for all transport managers.
 *
 *  The purpose of an untemplated base mesh is to simplify
 *  initialization of libraries, etc.
 */
class SOLVERS_EXPORT TransportManager
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  Initializes external libraries
   *
   *  @param argc       command line count
   *  @param argv       command line values
   */
  TransportManager(int argc, char *argv[])
    : d_flag(true)
  {
    callow_initialize(argc, argv);
  }

  /// Default constructor
  TransportManager()
    : d_flag(false)
  {
    /* ... */
  }

  /// Destructor.  This initializes all external libraries.
  virtual ~TransportManager()
  {
    if (d_flag) callow_finalize();
  }

private:

  /// Flag for initialization.
  bool d_flag;

};


} // end namespace detran

#endif // TRANSPORTMANAGER_HH_ 

//---------------------------------------------------------------------------//
//              end of file TransportManager.hh
//---------------------------------------------------------------------------//
