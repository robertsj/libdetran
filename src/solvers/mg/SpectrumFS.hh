//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumFS.hh
 *  @brief SpectrumFS class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SPECTRUMFS_HH_
#define detran_SPECTRUMFS_HH_

#include "SpectrumBase.hh"

namespace detran
{

/**
 *  @class SpectrumFS
 *  @brief Use the solution to a fixed source in an infinite homogeneous medium
 *
 *  We solve the material-dependent fixed source problem
 *  @f[
 *     (\mathbf{T} - \mathbf{S} - \frac{1}{k}\mathbf{XF}^T) \xi = b \, ,
 *
 *  @f]
 *  where @f$ b @f$ is a user-defined spectral forcing function that
 *  defaults to uniform across groups.
 */
class SpectrumFS: public SpectrumBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//


  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  SpectrumFS(SP_input input, SP_material material, SP_mesh mesh);

  virtual ~SpectrumFS(){}

  /// Produce the edit region-dependent spectra
  virtual vec2_dbl spectrum(double keff = 1.0);

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

};

} // end namespace detran

#endif /* detran_SPECTRUMFS_HH_ */

//----------------------------------------------------------------------------//
//              end of file SpectrumFS.hh
//----------------------------------------------------------------------------//
