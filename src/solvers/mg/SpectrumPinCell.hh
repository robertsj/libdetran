//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SpectrumPinCell.hh
 *  @brief SpectrumPinCell class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SPECTRUMPINCELL_HH_
#define detran_SPECTRUMPINCELL_HH_

#include "SpectrumBase.hh"
#include "geometry/PinCell.hh"

namespace detran
{

/**
 *  @class SpectrumPinCell
 *  @brief Use the critical spectrum of a pin cell for a region
 *
 *  This sceme allows a user to define representative 2-D pin cells
 *  for which to compute the critical eigenspectrum.  This spectrum is
 *  then applied to any fine mesh assigned the region for which
 *  the pin cell is defined.
 *
 *  The pin cell definitions are given in a separate pin cell
 *  input database named "mgpc_pincell_db".
 *
 */
class SpectrumPinCell: public SpectrumBase
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef std::vector<detran_geometry::PinCell::SP_pincell>   vec_pincell;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  SpectrumPinCell(SP_input input, SP_material material, SP_mesh mesh);

  virtual ~SpectrumPinCell(){}

  /// Produce the edit region-dependent spectra
  virtual vec2_dbl spectrum(double keff = 1.0);

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_input d_pincell_db;
  vec_pincell d_pincells;

};

} // end namespace detran

#endif /* detran_SPECTRUMPINCELL_HH_ */

//----------------------------------------------------------------------------//
//              end of file SpectrumPinCell.hh
//----------------------------------------------------------------------------//
