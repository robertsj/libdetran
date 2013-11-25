//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumGS.hh
 *  @brief SpectrumGS class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SPECTRUMGS_HH_
#define detran_SPECTRUMGS_HH_

#include "SpectrumBase.hh"

namespace detran
{

/**
 *  @class SpectrumGS
 *  @brief Use the dominant mode from the G-S eigenvalue problem as spectrum
 *
 *  Adams and Morel developed a two energy grid method for accelerating
 *  Gauss-Seidell upscatter interations.  The basic idea is to solve a
 *  one group diffusion problem for the error, using a special energy shape
 *  function to map between the one group and multigroup specta.  The shape
 *  they use is the eigenvector corresponding to the spectral radius of
 *  the Gauss-Seidel method for a homogeneous infinite medium, defined
 *  by
 *  @f[
 *     (\mathbf{T} - \mathbf{S}_L - \mathbf{S}_D)^{-1}
 *       \mathbf{S}_U \xi = \rho \xi \, .
 *
 *  @f]
 *  Here, @f$ \xi @f$ is the shape function of interest and is computed for
 *  each group.  Whereas Adams and Morel (and later, Evans, Clarno, and
 *  Morel) collapse to just one group, here, we use the shape function to
 *  map between the fine mesh and any allowable coarse energy grid.
 */
class SpectrumGS: public SpectrumBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//


  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  SpectrumGS(SP_input input, SP_material material, SP_mesh mesh);

  virtual ~SpectrumGS(){}

  /// Produce the edit region-dependent spectra
  virtual vec2_dbl spectrum(double keff = 1.0);

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//


};

} // end namespace detran

#endif /* detran_SPECTRUMGS_HH_ */

//----------------------------------------------------------------------------//
//              end of file SpectrumGS.hh
//----------------------------------------------------------------------------//
