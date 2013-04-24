//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SyntheticDiscreteSource.hh
 *  @brief  SyntheticDiscreteSource
 *  @author Jeremy Roberts
 *  @date   Nov 15, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_SYNTHETICDISCRETESOURCE_HH_
#define detran_SYNTHETICDISCRETESOURCE_HH_

#include "SyntheticSource.hh"

namespace detran
{

class KINETICS_EXPORT SyntheticDiscreteSource: public SyntheticSource
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef SyntheticSource       Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Mesh
   *  @param quadrature     Quadrature
   *  @param material       Material
   */
  SyntheticDiscreteSource(const size_t number_groups,
                          SP_mesh mesh,
                          SP_quadrature quadrature,
                          SP_material material);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SYNTHETIC SOURCES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Get moments source for cell.
   *
   *  Units are n/cc-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   */
  double source(const size_t cell,
                const size_t group);

  /**
   *  @brief Get discrete source for cell and cardinal angle.
   *
   *  Units are n/cc-ster-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   *  @param angle  Cardinal angle index
   */
  double source(const size_t cell,
                const size_t group,
                const size_t angle);

  /**
   *  @brief Build the synthetic source given previous iterates
   *
   *  Depending on the order of the method, we may use
   */
  virtual void build(const double dt,
                     const vec_states &states,
                     const vec_precursors &precursors,
                     const size_t order);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Discrete source [group][angle][space]
  vec3_dbl d_source;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "SyntheticDiscreteSource.i.hh"

#endif // detran_SYNTHETICDISCRETESOURCE_HH_

//---------------------------------------------------------------------------//
//              end of file SyntheticDiscreteSource.hh
//---------------------------------------------------------------------------//
