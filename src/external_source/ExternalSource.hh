//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ExternalSource.hh
 *  @brief ExternalSource class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_external_source_EXTERNALSOURCE_HH_
#define detran_external_source_EXTERNALSOURCE_HH_

#include "external_source/external_source_export.hh"
#include "angle/Quadrature.hh"
#include "geometry/Mesh.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_external_source
{

//----------------------------------------------------------------------------//
/**
 *  @class ExternalSource
 *  @brief Base volume source class
 *
 *  Consider the general transport equation in operator form:
 *  @f[
 *      \mathcal{L}_g \psi(\vec{r}, \hat{\Omega}, E_g)
 *        = Q(\vec{r}, \hat{\Omega}, E_g) \, .
 *  @f]
 *  The source \f$ Q \f$ is a function of space, angle, and energy.  Because
 *  sources will be used in different ways throughout Detran, a useful
 *  interface allows a client to get discrete angular sources (for use in
 *  SN or MOC calculations) and moment sources (for use in diffusion).  For
 *  the latter, only the isotropic component is available.
 *
 *  @note Source representation in moment form is limited to isotropic
 */
//----------------------------------------------------------------------------//

class EXTERNAL_SOURCE_EXPORT ExternalSource
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<ExternalSource>      SP_externalsource;
  typedef std::vector<SP_externalsource>            vec_externalsource;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef detran_utilities::size_t                  size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_groups  Number of energy groups
   *  @param mesh           Pointer to mesh
   *  @param quadrature     Pointer to angular quadrature
   *  @param discrete       Flag for discrete sources
   */
  ExternalSource(size_t         number_groups,
                 SP_mesh        mesh,
                 SP_quadrature  quadrature,
                 bool           discrete = false);

  /// Virtual destructor
  virtual ~ExternalSource(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  size_t number_groups() const { return d_number_groups; }
  bool   is_discrete() const { return d_discrete; }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EXTERNAL SOURCES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Get moments source for cell.
   *
   *  Units are n/cc-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   */
  virtual double source(const size_t cell,
                        const size_t group) = 0;

  /**
   *  @brief Get discrete source for cell and cardinal angle.
   *
   *  Units are n/cc-ster-sec
   *
   *  @param cell   Mesh cell
   *  @param group  Energy group
   *  @param angle  Cardinal angle index
   */
  virtual double source(const size_t cell,
                        const size_t group,
                        const size_t angle) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Cartesian mesh.
  SP_mesh d_mesh;
  /// Quadrature
  SP_quadrature d_quadrature;
  /// Number of groups.
  const size_t d_number_groups;
  /// Number of angles
  size_t d_number_angles;
  /// Am I ready?
  //int d_initialized;
  /// Discrete flag
  bool d_discrete;

};

EXTERNAL_SOURCE_TEMPLATE_EXPORT(detran_utilities::SP<ExternalSource>)
EXTERNAL_SOURCE_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<ExternalSource> >)

} // end namespace detran_external_source

#endif /* detran_external_sourceEXTERNALSOURCE_HH_ */

//----------------------------------------------------------------------------//
//              end of ExternalSource.hh
//----------------------------------------------------------------------------//
