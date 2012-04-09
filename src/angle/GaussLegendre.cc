//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussLegendre.cc
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Angle headers
#include "GaussLegendre.hh"

// System headers
#include <iostream>
#include <cstdio>
#include <cmath>

namespace detran
{

GaussLegendre::GaussLegendre(int order)
  : Quadrature(order,
               1,
               order,
               "GaussLegendre")
{
  Require(order % 2 == 0); // need an even order

  // Set GL values for one half-space based on the SN order specified.
  // For 2:2:20, we use hard-coded, high precision values, though I honestly
  // doubt the value of storing more than 16 decimal places.
  // For order > 20, the quadrature is computed numerically.
  switch (order)
  {
  case 2:
  {
    // GL S-2 Quadrature Set.
    d_mu[0] = 0.5773502691896257645091488;
    d_weight[0] = 1.0000000000000000000000000;
    break;
  }

  case 4:
  {
    // GL S-4 Quadrature Set
    d_mu[0] = 0.8611363115940525752239465;
    d_mu[1] = 0.3399810435848562648026658;
    d_weight[0] = 0.3478548451374538573730639;
    d_weight[1] = 0.6521451548625461426269361;
    break;
  }

  case 6:
  {
    // GL S-6 Quadrature Set
    d_mu[0] = 0.9324695142031520278123016;
    d_mu[1] = 0.6612093864662645136613996;
    d_mu[2] = 0.2386191860831969086305017;
    d_weight[0] = 0.1713244923791703450402961;
    d_weight[1] = 0.3607615730481386075698335;
    d_weight[2] = 0.4679139345726910473898703;
    break;
  }

  case 8:
  {
    // GL S-8 Quadrature Set
    d_mu[0] = 0.9602898564975362316835609;
    d_mu[1] = 0.7966664774136267395915539;
    d_mu[2] = 0.5255324099163289858177390;
    d_mu[3] = 0.1834346424956498049394761;
    d_weight[0] = 0.1012285362903762591525314;
    d_weight[1] = 0.2223810344533744705443560;
    d_weight[2] = 0.3137066458778872873379622;
    d_weight[3] = 0.3626837833783619829651504;
    break;
  }

  case 10:
  {
    // GL S-10 Quadrature Set
    d_mu[0] = 0.9739065285171717200779640;
    d_mu[1] = 0.8650633666889845107320967;
    d_mu[2] = 0.6794095682990244062343274;
    d_mu[3] = 0.4333953941292471907992659;
    d_mu[4] = 0.1488743389816312108848260;
    d_weight[0] = 0.0666713443086881375935688;
    d_weight[1] = 0.1494513491505805931457763;
    d_weight[2] = 0.2190863625159820439955349;
    d_weight[3] = 0.2692667193099963550912269;
    d_weight[4] = 0.2955242247147528701738930;
    break;
  }

  case 12:
  {
    // GL S-12 Quadrature Set */
    d_mu[0] = 0.9815606342467192506905491;
    d_mu[1] = 0.9041172563704748566784659;
    d_mu[2] = 0.7699026741943046870368938;
    d_mu[3] = 0.5873179542866174472967024;
    d_mu[4] = 0.3678314989981801937526915;
    d_mu[5] = 0.1252334085114689154724414;
    d_weight[0] = 0.0471753363865118271946160;
    d_weight[1] = 0.1069393259953184309602547;
    d_weight[2] = 0.1600783285433462263346525;
    d_weight[3] = 0.2031674267230659217490645;
    d_weight[4] = 0.2334925365383548087608499;
    d_weight[5] = 0.2491470458134027850005624;
    break;
  }

  case 14:
  {
    // GL S-14 Quadrature Set
    d_mu[0] = 0.9862838086968123388415973;
    d_mu[1] = 0.9284348836635735173363911;
    d_mu[2] = 0.8272013150697649931897947;
    d_mu[3] = 0.6872929048116854701480198;
    d_mu[4] = 0.5152486363581540919652907;
    d_mu[5] = 0.3191123689278897604356718;
    d_mu[6] = 0.1080549487073436620662447;
    d_weight[0] = 0.0351194603317518630318329;
    d_weight[1] = 0.0801580871597602098056333;
    d_weight[2] = 0.1215185706879031846894148;
    d_weight[3] = 0.1572031671581935345696019;
    d_weight[4] = 0.1855383974779378137417166;
    d_weight[5] = 0.2051984637212956039659241;
    d_weight[6] = 0.2152638534631577901958764;
    break;
  }

  case 16:
  {
    // GL S-16 Quadrature Set
    d_mu[0] = 0.9894009349916499325961542;
    d_mu[1] = 0.9445750230732325760779884;
    d_mu[2] = 0.8656312023878317438804679;
    d_mu[3] = 0.7554044083550030338951012;
    d_mu[4] = 0.6178762444026437484466718;
    d_mu[5] = 0.4580167776572273863424194;
    d_mu[6] = 0.2816035507792589132304605;
    d_mu[7] = 0.0950125098376374401853193;
    d_weight[0] = 0.0271524594117540948517806;
    d_weight[1] = 0.0622535239386478928628438;
    d_weight[2] = 0.0951585116824927848099251;
    d_weight[3] = 0.1246289712555338720524763;
    d_weight[4] = 0.1495959888165767320815017;
    d_weight[5] = 0.1691565193950025381893121;
    d_weight[6] = 0.1826034150449235888667637;
    d_weight[7] = 0.1894506104550684962853967;
    break;
  }

  case 18:
  {
    // GL S-18 Quadrature Set
    d_mu[0] = 0.9915651684209309467300160;
    d_mu[1] = 0.9558239495713977551811959;
    d_mu[2] = 0.8926024664975557392060606;
    d_mu[3] = 0.8037049589725231156824175;
    d_mu[4] = 0.6916870430603532078748911;
    d_mu[5] = 0.5597708310739475346078715;
    d_mu[6] = 0.4117511614628426460359318;
    d_mu[7] = 0.2518862256915055095889729;
    d_mu[8] = 0.0847750130417353012422619;
    d_weight[0] = 0.0216160135264833103133427;
    d_weight[1] = 0.0497145488949697964533349;
    d_weight[2] = 0.0764257302548890565291297;
    d_weight[3] = 0.1009420441062871655628140;
    d_weight[4] = 0.1225552067114784601845191;
    d_weight[5] = 0.1406429146706506512047313;
    d_weight[6] = 0.1546846751262652449254180;
    d_weight[7] = 0.1642764837458327229860538;
    d_weight[8] = 0.1691423829631435918406565;
    break;
  }

  case 20:
  {
    // GL S-20 Quadrature Set
    d_mu[0] = 0.9931285991850949247861224;
    d_mu[1] = 0.9639719272779137912676661;
    d_mu[2] = 0.9122344282513259058677524;
    d_mu[3] = 0.8391169718222188233945291;
    d_mu[4] = 0.7463319064601507926143051;
    d_mu[5] = 0.6360536807265150254528367;
    d_mu[6] = 0.5108670019508270980043641;
    d_mu[7] = 0.3737060887154195606725482;
    d_mu[8] = 0.2277858511416450780804962;
    d_mu[9] = 0.0765265211334973337546404;
    d_weight[0] = 0.0176140071391521183118620;
    d_weight[1] = 0.0406014298003869413310400;
    d_weight[2] = 0.0626720483341090635695065;
    d_weight[3] = 0.0832767415767047487247581;
    d_weight[4] = 0.1019301198172404350367501;
    d_weight[5] = 0.1181945319615184173123774;
    d_weight[6] = 0.1316886384491766268984945;
    d_weight[7] = 0.1420961093183820513292983;
    d_weight[8] = 0.1491729864726037467878287;
    d_weight[9] = 0.1527533871307258506980843;
    break;
  }

  default:
  {
    // generate the parameters on-the-fly
    double tmp_mu[d_number_angles_octant];
    double tmp_wt[d_number_angles_octant];
    generate_parameters(order, tmp_mu, tmp_wt);
    for (int i = 0; i < d_number_angles_octant; i++)
    {
      d_mu[i]     = tmp_mu[i];
      d_weight[i] = tmp_wt[i];
    }
    break;
  }
  }

} // end constructor

/*!
 * \brief Generate Gauss-Legendre parameters.
 *
 * Modified version of the function gauleg from <B> Numerical Recipes in C</B>
 * (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
 * & B.P. Flannery
 *
 * \param order     number of points (i.e. the quadrature order)
 * \param mu        temporary array for Gauss points
 * \param wt        temporary array for weights
 *
 */
void GaussLegendre::generate_parameters(int order, double *mu, double *wt)
{
  int m, j, i;
  double z1, z, xm, xl, pp, p3, p2, p1;
  // The roots are symmetric, so we only find half of them.
  m = (order + 1) / 2;
  for (i = 1; i <= m; i++)
  { /* Loop over the desired roots. */
    z = cos(pi * (i - 0.25) / (order + 0.5));
    // Starting with the above approximation to the ith root, we enter
    // the main loop of refinement by Newton's method.
    int count = 0;
    do
    {
      p1 = 1.0;
      p2 = 0.0;
      // Recurrence to get Legendre polynomial.
      for (j = 1; j <= order; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      // p1 is now the desired Legendre polynomial. We next compute
      // pp, its derivative, by a standard relation involving also
      // p2, the polynomial of one lower order.
      pp = order * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp; // <-- Newton's method
    } while (std::abs(z - z1) > 1e-15 && ++count < 200);
    // Store the root and compute the weight.
    mu[i - 1] = z;
    wt[i - 1] = 2.0 / ((1.0 - z * z) * pp * pp);
  }
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Display the quadrature.
 */
void GaussLegendre::display() const
{

  using std::cout;
  using std::endl;
  using std::printf;


  cout << endl;
  cout << d_name << " abscissa and weights: " << endl << endl;
  cout << "   m         mu              wt       " << endl;
  cout << "  ---   --------------  --------------" << endl;

  double weight_sum = 0;
  for (int i = 0; i < d_number_angles_octant; ++i)
  {
    printf("%4i    %13.10f   %13.10f \n", i, d_mu[i], d_weight[i]);
    weight_sum += d_weight[i];
  }
  cout << endl << "  The sum of the weights is " << weight_sum << endl;
  cout << endl;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//                 end of GaussLegendre.cc
//---------------------------------------------------------------------------//
