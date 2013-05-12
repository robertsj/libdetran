//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Material.hh
 *  @author Jeremy Roberts
 *  @brief  Material class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_material_MATERIAL_HH_
#define detran_material_MATERIAL_HH_

#include "material/material_export.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <string>
#ifdef DETRAN_ENABLE_BOOST
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif

namespace detran_material
{

//---------------------------------------------------------------------------//
/**
 *  @class Material
 *  @brief Simple cross section container.
 *
 *  All data is stored with the material index changing fastest.  This
 *  appears to be the best storage scheme with respect to memory access.
 */
//---------------------------------------------------------------------------//
class MATERIAL_EXPORT Material
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Material> SP_material;
  typedef detran_utilities::vec_dbl      vec_dbl;
  typedef detran_utilities::vec2_dbl     vec2_dbl;
  typedef detran_utilities::vec3_dbl     vec3_dbl;
  typedef detran_utilities::vec_int      vec_int;
  typedef detran_utilities::vec2_int     vec2_int;
  typedef detran_utilities::vec_size_t   vec_size_t;
  typedef detran_utilities::vec2_size_t  vec2_size_t;
  typedef detran_utilities::size_t       size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param    number_materials    Number of materials.
   *  @param    number_groups       Number of energy groups.
   *  @param    downscatter         Switch on to use only downscatter.
   */
  Material(const size_t number_materials,
           const size_t number_groups,
           std::string  name = "no name given");

  /// Virtual destructor
  virtual ~Material(){}

  /// SP constructor
  static SP_material Create(const size_t number_materials,
                            const size_t number_groups,
                            std::string  name = "no name given");

  //--------------------------------------------------------------------------//
  // Setters
  //--------------------------------------------------------------------------//

  /**
   *  @brief Explicitly turn on downscatter-only
   */
  void set_downscatter(bool v, bool tran = false)
  {
    if (tran) d_downscatter[1] = v;
    d_downscatter[0] = v;
    if (d_finalized) finalize();
  }

  void set_sigma_t(size_t m, size_t g, double v);
  void set_sigma_a(size_t m, size_t g, double v);
  void set_nu_sigma_f(size_t m, size_t g, double v);
  void set_sigma_f(size_t m, size_t g, double v);
  void set_nu(size_t m, size_t g, double v);
  void set_chi(size_t m, size_t g, double v);
  void set_sigma_s(size_t m, size_t g, size_t gp, double v);
  void set_diff_coef(size_t m, size_t g, double v);

  // Vectorized setters

  void set_sigma_t(size_t m, vec_dbl &v);
  void set_sigma_a(size_t m, vec_dbl &v);
  void set_nu_sigma_f(size_t m, vec_dbl &v);
  void set_sigma_f(size_t m, vec_dbl &v);
  void set_nu(size_t m, vec_dbl &v);
  void set_chi(size_t m, vec_dbl &v);
  void set_sigma_s(size_t m, size_t g, vec_dbl &v);
  void set_diff_coef(size_t m, vec_dbl &v);

  //------------------------------------------------------------------------//
  // Getters
  //------------------------------------------------------------------------//

  virtual double sigma_t(size_t m, size_t g) const;
  virtual double sigma_a(size_t m, size_t g) const;
  virtual double nu_sigma_f(size_t m, size_t g) const;
  virtual double sigma_f(size_t m, size_t g) const;
  virtual double nu(size_t m, size_t g) const;
  virtual double chi(size_t m, size_t g) const;
  virtual double sigma_s(size_t m, size_t g, size_t gp) const;
  virtual double diff_coef(size_t m, size_t g) const;

  // Vectorized getters

  virtual vec_dbl sigma_t(size_t m) const;
  virtual vec_dbl sigma_a(size_t m) const;
  virtual vec_dbl nu_sigma_f(size_t m) const;
  virtual vec_dbl sigma_f(size_t m) const;
  virtual vec_dbl nu(size_t m) const;
  virtual vec_dbl chi(size_t m) const;
  virtual vec2_dbl sigma_s(size_t m) const;
  virtual vec_dbl diff_coef(size_t m) const;

  //------------------------------------------------------------------------//
  // OTHER ACCESSORS
  //------------------------------------------------------------------------//

  size_t number_groups() const
  {
    return d_number_groups;
  }

  size_t number_materials() const
  {
    return d_number_materials;
  }

  /**
   *  @brief Lower scatter group bound.
   *
   *  This is the *lowest* index (highest energy) \f$ g' \f$
   *  that leads to downscatter for a given outgoing group \f$ g \f$.
   *
   *  @param g        Row of the scattering matrix
   *  @param tran     Flag for accessing transpose of S
   */
  size_t lower(size_t g, bool tran = false) const;

  /**
   *  @brief Upper scatter group bound.
   *
   *  This is the *highest* index (lowest energy) \f$ g' \f$
   *  that upscatters size_to the outgoing group \f$ g \f$.
   *
   *  @param g        Row of the scattering matrix
   *  @param tran     Flag for accessing transpose of S
   */
  size_t upper(size_t g, bool tran = false) const;

  /// Do we do only downscatter?
  bool downscatter(bool tran = false) const;

  /**
   *  @brief Index below which upscatter doesn't occur for any material.
   *
   *  For adjoint problems, this is the group above which
   */
  size_t upscatter_cutoff(bool tran = false) const;

  /**
   *  @brief Compute the absorption cross section from total and scattering.
   *  @note this overwrites any data for \f$ \Sigma_a \f$ already stored.
   */
  void compute_sigma_a();

  /**
   *  @brief Compute the diffusion coefficient from \f$ \Sigma_t \f$.
   *
   *  Assuming isotropic scattering in the LAB, the diffusion
   *  coefficient is simply \f$ D = 1/3\Sigma_t \f$.
   *
   *  @todo Update diffusion definition if anisotropic scattering
   *        is added.
   *  @note This overwrites any data for \f$ D \f$ already stored.
   */
  void compute_diff_coef();

  /// Computes scattering bounds and absorption cross section.
  void finalize();

  /// Pretty print the material database.
  virtual void display();

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Material name
  std::string d_name;
  /// Number of groups
  size_t d_number_groups;
  /// Number of materials
  size_t d_number_materials;
  /// Downscatter switch (when true, upscatter ignored)
  bool d_downscatter[2];
  /// Total cross section [material, group]
  vec2_dbl d_sigma_t;
  /// Absorption cross section [material, group]
  vec2_dbl d_sigma_a;
  /// nu * Fission [material, group]
  vec2_dbl d_nu_sigma_f;
  /// Fission [material, group]
  vec2_dbl d_sigma_f;
  /// nu [material, group]
  vec2_dbl d_nu;
  /// Fission spectrum [material, group]
  vec2_dbl d_chi;
  /// Scatter [material, group<-, group']
  vec3_dbl d_sigma_s;
  /// Diffusion coefficient [material, group]
  vec2_dbl d_diff_coef;
  /// Scatter bounds applied to all materials [group, 2]
  vec2_size_t d_scatter_bounds;
  /// Groups equal to or above cutoff are subject to upscatter iterations
  size_t d_upscatter_cutoff[2];
  /// Are we ready to be used?
  bool d_finalized;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  void material_display();

#ifdef DETRAN_ENABLE_BOOST

  /// Default constructor needed for serialization
  Material(){}

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_number_groups;
    ar & d_number_materials;
    ar & d_downscatter;
    ar & d_sigma_t;
    ar & d_sigma_a;
    ar & d_nu_sigma_f;
    ar & d_sigma_f;
    ar & d_nu;
    ar & d_chi;
    ar & d_sigma_s;
    ar & d_diff_coef;
    ar & d_scatter_bounds;
    ar & d_upscatter_cutoff;
    ar & d_finalized;
  }

#endif

};

MATERIAL_TEMPLATE_EXPORT(detran_utilities::SP<Material>)

} // end namespace detran_material

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Material.i.hh"

#endif /* detran_material_MATERIAL_HH_ */
