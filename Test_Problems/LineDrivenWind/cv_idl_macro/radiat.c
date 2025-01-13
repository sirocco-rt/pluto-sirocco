/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Computes the right-hand side for tabulated cooling.

  \author Amin Mosallanezhad
  \date   December 31, 2024
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"

extern double gCooling_x1, gCooling_x2, gCooling_x3;
extern double gSirocco_xi, gSirocco_t_r;
extern double gComp_heat_pre, gComp_cool_pre;
extern double gXray_heat_pre, gLine_cool_pre;
extern double gBrem_cool_pre, gXi_ion_pre;

/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 *
 ******************************************************************* */
{

  double T, prs, mu, lx, r, xi, tx;
  double rho, nH, ne;
  static double E_cost;

  double sqT, sqxi, sqsqxi;

  double comp_heat, comp_cool, xray_heat;
  double line_cool, brem_cool;
  double lambda;                     // Net heating/cooling rate
  double ne_ratio;

  g_minCoolingTemp = 1.e4;           // Minimum cooling temperature (in Kelvin)

  mu = g_inputParam[MU];                           // Mean particle mass
  lx = g_inputParam[L_star] * g_inputParam[f_x];   // X-ray luminosity
  r  = gCooling_x1 * UNIT_LENGTH;                  // Radius in real units


  rho = v[RHO] * UNIT_DENSITY;       // Density in physical units

  /* Compute hydrogen number density assuming stellar abundances */
  nH = rho / (1.43 * CONST_mp);  // 1.43 * proton mass for stellar composition
  // ne = 1.21 * nH;

  E_cost = UNIT_LENGTH / UNIT_DENSITY / pow(UNIT_VELOCITY, 3.0);

/* ---------------------------------------------
            Get pressure and temperature
   --------------------------------------------- */

  prs = v[RHOE] * (g_gamma - 1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs / (g_gamma - 1.0);
  }


  T  = prs / v[RHO] * KELVIN * mu;

  if (T != T){
    printf ("! Radiat(): Nan found: rho = %12.6e, prs = %12.6e\n", v[RHO], prs);
    QUIT_PLUTO(1);
  }


  if (T < g_minCoolingTemp) {
    rhs[RHOE] = 0.0;
    return;
  }



  /* Initialize ionization parameter xi if necessary */
  if (g_time < 3.0) {
    xi = lx / nH / r / r;  // Initial ionization parameter
    tx = g_inputParam[T_x];
  } else {
    // xi = sirocco_xi[k][j][i] * xi_ion_pre;  // Update ionization parameter
    xi = gSirocco_xi;                          // Update ionization parameter
    tx = gSirocco_t_r;
  }


  sqT    = sqrt(T);           // T^(0.5)
  sqxi   = sqrt(xi);          // xi^(0.5)
  sqsqxi = pow(xi, 0.25);      // xi^(0.25)


  /* Compute electron number density using ne_rat function */
  #if PY_CONNECT
    ne = nH * ne_rat(T);
  #else
    /* Assume fully ionized gas */
    ne = 1.21 * nH;
  #endif


// Pure Blondin rates without prefactors!

  // comp_heat = (8.9e-36 * xi * tx );
  // comp_cool = (8.9e-36 * xi * (4.0 * T));
  // xray_heat = (1.5e-21 * (sqsqxi / sqT) * (1 - (T / tx)));
  // line_cool = ((1e-16 * exp(-1.3e5 / T) / sqxi / T) +
  //           fmin(fmin(1e-24, 5e-27 * sqT), 1.5e-17 / T));
  // brem_cool = (3.3e-27 * sqT);



// Blondin rates with prefactors!

  comp_heat = gComp_heat_pre * (8.9e-36 * xi * tx);
  comp_cool = gComp_cool_pre * (8.9e-36 * xi * (4.0 * T));
  xray_heat = gXray_heat_pre * (1.5e-21 * (sqsqxi / sqT));
  line_cool = gLine_cool_pre * ((1e-16 * exp(-1.3e5 / T) / sqxi / T) +
              fmin(fmin(1e-24, 5e-27 * sqT), 1.5e-17 / T));
  brem_cool = gBrem_cool_pre * (3.3e-27 * sqT);


  lambda = nH * (ne * comp_heat + nH * xray_heat
              -  ne * comp_cool - ne * line_cool - ne * brem_cool);

  if (v[TRC] > 0.01){
    rhs[RHOE] = lambda * E_cost;
  } else {
    rhs[RHOE] = 0.0;
  }


/* ----------------------------------------------
    Temperature cutoff
   ---------------------------------------------- */

  rhs[RHOE] *= 1.0 - 1.0 /cosh( pow( T / g_minCoolingTemp, 12));

}
