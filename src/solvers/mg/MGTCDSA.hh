//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGTCDSA.hh
 *  @brief MGTCDSA class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGTCDSA_HH_
#define detran_MGTCDSA_HH_

#include "MGTransportOperator.hh"
#include "MGDSA.hh"

namespace detran
{

/**
 *  @class MGTCDSA
 *  @brief Multigroup transport-corrected diffusion preconditioner
 *
 *  Starting from a diffusion-based preconditioner, the Newton-Shulz
 *  method is used to apply a correction based on zero or more
 *  applications of an approximate transport operator in conjunction with
 *  application of the full transport operator.  Two or three approximate
 *  transport applications seems to yield a good start for the final
 *  full transport correction.
 */
template <class D>
class MGTCDSA: public MGPreconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef Base::SP_preconditioner           SP_base;
  typedef detran_utilities::SP<MGTCDSA>     SP_pc;
  typedef DiffusionLossOperator             DSA_T;
  typedef MGTransportOperator<D>            MGTO_T;
  typedef typename MGTO_T::SP_operator      SP_mgto;
  typedef MGDSA::SP_preconditioner          SP_mgdsa;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input            Input database
   *  @param material         Material database
   *  @param mesh             Cartesian mesh
   *  @param ssource          Scattering source
   *  @param fsource          Fission source
   *  @param P                Original preconditioner
   *  @param A                Full transport operator
   *  @param cutoff           First group included in solve
   *  @param include_fission  Flag to treat fission like scatter
   *  @param adjoint          Flag for adjoint problems
   */
  MGTCDSA(SP_input                input,
          SP_material             material,
          SP_mesh                 mesh,
          SP_scattersource        ssource,
          SP_fissionsource        fsource,
          SP_base                 P,
          SP_mgto                 A,
          size_t                  cutoff,
          bool                    include_fission,
          bool                    adjoint);

  /// virtual destructor
  virtual ~MGTCDSA(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// Solve Px = b
  void apply(Vector &b, Vector &x);

  /// build the preconditioner for a new keff
  void build(const double keff = 1.0, SP_state state = SP_state(0));

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Full transport operator
  SP_mgto d_A;
  /// Approximate transport operator
  SP_mgto d_A_tilde;
  /// Original (diffusion) preconditioner
  SP_base d_P;
  /// Number of low-order corrections
  size_t d_number_coarse_corrections;
  /// Use the high order correction
  bool d_use_fine_correction;
  /// Count of A applications
  size_t d_A_count;
  /// Count of A_tilde applications
  size_t d_A_tilde_count;
  /// Count of P applications
  size_t d_P_count;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  void apply(Vector &b, Vector &x, const size_t k);

};

} // end namespace detran

#endif /* detran_MGTCDSA_HH_ */

//----------------------------------------------------------------------------//
//              end of file MGTCDSA.hh
//----------------------------------------------------------------------------//
