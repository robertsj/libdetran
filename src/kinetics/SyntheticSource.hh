//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SyntheticSource.hh
 *  @brief  SyntheticSource
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_SYNTHETICSOURCE_HH_
#define detran_SYNTHETICSOURCE_HH_

#include "KineticsMaterial.hh"
#include "Precursors.hh"
#include "external_source/ExternalSource.hh"
#include "transport/State.hh"
#include "utilities/DBC.hh"
#include "BDFCoefficients.hh"

namespace detran
{

/**
 *  @class SyntheticSource
 *  @brief Time step source contribution from previous step solution
 *
 *  In all the time stepping schemes implemented in Detran,
 *  the time step is cast as a fixed source problem with
 *  a synthetic source.  This class builds the source from
 *  previous state vectors (using either the angular or scalar
 *  flux) and previous precursor vectors.
 *
 *  This class implements synthetic sources based on
 *  BDF discretizations of orders 1-6.  Note, first order
 *  BDF is equivalent to backward Euler and is used (via a half-step)
 *  to implement the implicit midpoint rule as used in
 *  PARTISN.
 *
 */
class KINETICS_EXPORT SyntheticSource:
  public detran_external_source::ExternalSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_external_source::ExternalSource  Base;
  typedef Base::SP_externalsource                 SP_base;
  typedef detran_utilities::SP<SyntheticSource>   SP_externalsource;
  typedef KineticsMaterial::SP_material           SP_material;
  typedef State::SP_state                         SP_state;
  typedef std::vector<SP_state>                   vec_states;
  typedef Precursors::SP_precursors               SP_precursors;
  typedef std::vector<SP_precursors>              vec_precursors;
  typedef detran_utilities::vec_dbl               vec_dbl;
  typedef detran_utilities::vec2_dbl              vec2_dbl;
  typedef detran_utilities::vec3_dbl              vec3_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Mesh
   *  @param quadrature     Angular mesh
   *  @param material       Material
   *  @param discrete       Flag for discrete source
   */
  SyntheticSource(const size_t  number_groups,
                  SP_mesh       mesh,
                  SP_quadrature quadrature,
                  SP_material   material,
                  bool          discrete = false);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SYNTHETIC SOURCES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Build the synthetic source given previous iterates
   *
   *  Depending on the order of the method, we may use
   */
  virtual void build(const double dt,
                     const vec_states &states,
                     const vec_precursors &precursors,
                     const size_t order = 1) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Material database
  SP_material d_material;

  /// Angular norm
  double d_norm;

};

KINETICS_TEMPLATE_EXPORT(detran_utilities::SP<SyntheticSource>)

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "SyntheticSource.i.hh"

#endif // SYNTHETICSOURCE_HH_ 

//---------------------------------------------------------------------------//
//              end of file SyntheticSource.hh
//---------------------------------------------------------------------------//
