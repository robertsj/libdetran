//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_MOC.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_SC_MOC inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_SC_MOC_I_HH_
#define detran_EQUATION_SC_MOC_I_HH_

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
inline void Equation_SC_MOC::solve(const size_t       region,
                                   const double       length,
                                   const double       width,
                                   moments_type      &source,
                                   double            &psi_in,
                                   double            &psi_out,
                                   moments_type      &phi,
                                   angular_flux_type &psi)
{
  using std::cout;
  using std::endl;

  // Preconditions.
  Require(region < d_mesh->number_cells());

  double sigma = d_material->sigma_t(d_mat_map[region], d_g);
  double length_over_sin = length * d_inv_sin[d_polar];
  double inv_volume = 1.0 / d_mesh->volume(region);

  // Coefficients from Hebert.
  double A = std::exp(-sigma * length_over_sin);
  double B = (1.0 - A) / sigma;
  double C = (length_over_sin / sigma) * (1.0 - B / length_over_sin);
  //double C = (length_over_sin / sigma) * (1.0 - (1.0-A)/(length_over_sin*sigma));

  // Segment outgoing angular flux.
  psi_out = A * psi_in + B * source[region];


//  if (!(source[region] > 0.0 ? psi_out > 0.0 : true))
//  {
//    cout << "  region = " << region << endl;
//    cout << "  psi_in = " << psi_in << endl;
//    cout << " psi_out = " << psi_out << endl;
//    cout << "   sigma = " << sigma << endl;
//    cout << " L / sin = " << length_over_sin << endl;
//    cout << "       A = " << A << endl;
//    cout << "       B = " << B << endl;
//    cout << "  source = " << source[region] << endl;
//  }
//  Ensure(source[region] > 0.0 ? psi_out > 0.0 : true);

  // Segment average angular flux.
  double psi_average = (1.0 / length_over_sin) * (B * psi_in + C * source[region]);

//  cout << "        t = " << length_over_sin << endl;
//  cout << "        A = " << A << endl;
//  cout << "        B = " << B << endl;
//  cout << "        C = " << C << endl;
//  cout << "   source = " << source[region] << endl;
//  cout << "  psi_in  = " << psi_in  << endl;
//  cout << "  psi_out = " << psi_out << endl;
//  cout << "  psi_seg = " << psi_average << endl;
//  cout << "   length = " << length << endl;
//  cout << "    space = " << d_spacing << endl;

  // Contribution to region average angular flux.
  psi_average *= width * length * inv_volume;

  // Contribution to region average scalar flux.
  phi[region] += d_quadrature->weight(d_angle) * psi_average;

//  cout << "  psi_reg = " << psi_average << endl;
//  cout << "   weight = " << d_quadrature->weight(d_angle) << endl;
//  cout << "  phi_reg = " << phi[region] << endl;

  // Store angular flux if needed.
  if (d_update_psi) psi[region] += psi_average;

}

} // end namespace detran

#endif /* detran_EQUATION_SC_MOC_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_MOC.i.hh
//---------------------------------------------------------------------------//
