//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   KineticsMaterial.hh
 *  @author robertsj
 *  @date   Oct 2, 2012
 *  @brief  KineticsMaterial class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_KINETICSMATERIAL_HH_
#define detran_KINETICSMATERIAL_HH_

#include "kinetics/kinetics_export.hh"
#include "material/Material.hh"

namespace detran
{

/**
 *  @class KineticsMaterial
 *  @brief Extends Material class for use in kinetics problems
 *
 *  The basic structure of the materials class for time-dependent
 *  problems is to expose the typical Material interface, extend
 *  it to include basic kinetics parameters, and to provide
 *  an interface for defining materials at a time.  For all
 *  cases, base cross section data is defined.  For a few
 *  parameters, a synthetic value is defined that accounts
 *  for the time discretization.  While in the most general
 *  case, materials will be dependent on the fine mesh, a
 *  region base material is likely present.  For the case of
 *  a fine mesh-dependent material, the user will have to
 *  explicitly create a material map with unique id's for each
 *  fine mesh.
 *
 *
 *  Time-dependent transport in reactor applications requires
 *  tracking of delayed neutrons.  This class contains the
 *  delayed neutron fractions (i.e. \f$ \beta_i \f$) for each
 *  delayed precursor group \f$ i \f$ for every material.
 *
 *  Some assumptions:
 *    - assume that \f$ v_g \f$ is time-independent
 *    - assume that \f$ \beta_i \f$ is material-dependent
 *    - assume that \f$ \lambda_i \f$ is material-independent
 *
 */
class KINETICS_EXPORT KineticsMaterial: public detran_material::Material
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<KineticsMaterial>    SP_material;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_materials           Number of materials
   *  @param number_energy_groups       Number of energy groups
   *  @param number_precursor_groups    Number of precursor groups
   *  @param name                       Material name
   */
  KineticsMaterial(const size_t number_materials,
                   const size_t number_energy_groups,
                   const size_t number_precursor_groups,
                   std::string  name = "KineticsMaterial");

  /// Virtual destructor
  virtual ~KineticsMaterial(){}

  /// SP constructor
  static SP_material Create(const size_t number_materials,
                            const size_t number_energy_groups,
                            const size_t number_precursor_groups,
                            std::string  name = "KineticsMaterial");

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Set the velocity for  group \f$ i \f$
  void set_velocity(const size_t g, const double v);
  /// Set the decay constant for group \f$ i \f$
  void set_lambda(const size_t i, const double v);
  /// Set delayed neutron fraction for a material and group
  void set_beta(const size_t m, const size_t i, const double v);
  /// Set delayed neutron fraction for a group (constant for all materials)
  void set_beta(const size_t i, const double v);
  /// Set delayed neutron spectrum
  void set_chi_d(const size_t m, const size_t i, const size_t g, double v);

  /// Get the velocity for group \f$ g \f$
  double velocity(const size_t g) const;
  /// Get the decay constant for group \f$ i \f$
  double lambda(const size_t i) const;
  /// Get delayed neutron fraction for a material and group
  virtual double beta(const size_t m, const size_t i) const;
  /// Get total delayed neutron fraction for a material
  virtual double beta_total(const size_t m) const;
  /// Get delayed neutron spectrum
  virtual double chi_d(const size_t m, const size_t i, const size_t g) const;

  /// Get the number of precursor groups
  size_t number_precursor_groups() const;

  /// Check that all values are positive, etc.
  void finalize();

  /// Display
  void display();

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Number of precursor groups
  size_t d_number_precursor_groups;
  /// Group velocities
  vec_dbl d_velocity;
  /// Decay constants [number precursors groups]
  vec_dbl d_lambda;
  /// Delayed neutron fractions [number precursor groups][number materials]
  vec2_dbl d_beta;
  /// Delayed neutron spectra [num e groups][num p groups][number materials]
  vec3_dbl d_chi_d;
  /// Verified
  bool d_verified;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  void kinetics_display();
};

KINETICS_TEMPLATE_EXPORT(detran_utilities::SP<KineticsMaterial>)
KINETICS_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<KineticsMaterial> >)

} // end namespace detran

#include "KineticsMaterial.i.hh"

#endif /* detran_KINETICSMATERIAL_HH_ */
