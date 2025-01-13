/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implements a source step using BLONDIN cooling.

  \details This file contains the implementation for advancing a
           source step in the simulation, incorporating the BLONDIN
           cooling mechanism to account for radiative cooling effects.

  \author Nick Higginbottom
  \date   July 28, 2018

  \revised_by A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date       December 31, 2024
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"

#define ITMAX 100
#define EPS 3.0e-8

/* Global variables */
double xi, sqxi, sqsqxi;       // Ionization parameters
double ne, nH, n;              // Number densities
double tx, E, hc_init, rho;    // Temperature, energy, initial heating/cooling rate, density
double dt_share;               // Shared time step in physical units
double comp_c_pre, comp_h_pre; // Precomputed heating/cooling coefficients
double line_c_pre, brem_c_pre, xray_h_pre, xi_ion_pre; // Precomputed coefficients
int flag;                      // General-purpose flag

/* Function prototypes */
double heatcool(double T);
double zfunc(double temp);
double zbrent(double (*func)(double), double x1, double x2, double tol);
double ne_rat(double temp);

/* ********************************************************************* */
/*!
 * \brief  Takes a source step using BLONDIN cooling.
 *
 * \param [in,out]  VV    Pointer to the PLUTO 3D data array containing primitive variables.
 * \param [in]      data  Pointer to the data structure containing additional variables.
 * \param [in]      dt    The current integration time step.
 * \param [in]      Dts   Pointer to the Time_Step structure.
 * \param [in]      grid  Pointer to an array of Grid structures.
 *
 *********************************************************************** */
void BlondinCooling(Data_Arr VV, const Data *data, double dt, timeStep *Dts, Grid *grid)
{
  int i, j, k;                       // Loop indices
  double p, T, p_f, T_f;             // Pressure and temperature variables
  double r;                          // Radius
  double t_l, t_u;                   // Lower and upper bounds for root finding
  double test;                       // Variable to test bracketing
  double E_f;                        // Final internal energy
  double hc_final;                   // Final heating/cooling rate
  double T_test;                     // Test temperature
  double mu;                         // Mean molecular weight
  double lx;                         // X-ray luminosity
  int reset_flag = 0;                // Reset flag

  dt_share = dt * UNIT_TIME;         // Time step in physical units for root finder
  lx = g_inputParam[L_star] * g_inputParam[f_x];  // X-ray luminosity
  mu = g_inputParam[MU];                          // Mean particle mass
  flag = 0;
  g_minCoolingTemp = 1.e4;           // Minimum cooling temperature (in Kelvin)

  /* Loop over the computational domain, excluding ghost zones */
  DOM_LOOP(k, j, i) {

    /* Retrieve precomputed heating and cooling coefficients from the data structure */
    comp_c_pre = data->comp_c_pre[k][j][i];
    comp_h_pre = data->comp_h_pre[k][j][i];
    line_c_pre = data->line_c_pre[k][j][i];
    brem_c_pre = data->brem_c_pre[k][j][i];
    xray_h_pre = data->xray_h_pre[k][j][i];
    xi_ion_pre = data->xi_ion_pre[k][j][i];

    /* Compute radius in physical units */
    r = grid->x[IDIR][i] * UNIT_LENGTH;  // Radius in real units

    /* Compute density and pressure  */
    rho = VV[RHO][k][j][i] * UNIT_DENSITY;  // Density in physical units
    p   = VV[PRS][k][j][i];                 // Pressure in code units

    /* Compute initial temperature in Kelvin */
    T = p / VV[RHO][k][j][i] * KELVIN * mu;

    /* Compute internal energy in physical units */
    E = (p * UNIT_PRESSURE) / (g_gamma - 1);

    /* Compute hydrogen number density assuming stellar abundances */
    nH = rho / (1.43 * CONST_mp);  // 1.43 * proton mass for stellar composition


    /* Initialize ionization parameter xi if necessary */
    if (g_time <= 3.0) {
      xi = lx / nH / r / r;  // Initial ionization parameter
      tx = g_inputParam[T_x];
    } else {
      // xi = sirocco_xi[k][j][i] * xi_ion_pre;  // Update ionization parameter
      xi = data->sirocco_xi[k][j][i];   // Update ionization parameter
      tx = data->sirocco_t_r[k][j][i];
    }


    /* Compute particle number density */
    n = rho / (mu * CONST_mp);

    /* Recompute temperature from energy */
    T = E * (2.0 / 3.0) / (n * CONST_kB);

    /* Skip if temperature is below minimum cooling temperature */
    if (T < g_minCoolingTemp) continue;

    /* Precompute xi powers for efficiency */
    sqxi   = sqrt(xi);       // xi^(0.5)
    sqsqxi = pow(xi, 0.25);

    /* Get initial heating/cooling rate */
    hc_init = heatcool(T);

    /* Bracket the solution temperature */
    t_l = T * 0.9;
    t_u = T * 1.1;
    test = zfunc(t_l) * zfunc(t_u);

    /* Adjust bounds until the root is bracketed */
    while (test > 0 && test == test) {
      t_l *= 0.9;
      t_u *= 1.1;
      test = zfunc(t_l) * zfunc(t_u);
    }

    if (test != test) {
      /* If test is NaN, cannot bracket the temperature */
      printf("Search has errored\n");
      T_f = T;
    } else {
      /* Root finding to solve for final temperature T_f */
      T_f = zbrent(zfunc, t_l, t_u, 1.0);
      hc_final = heatcool(T_f);

      /* Check if equilibrium temperature is crossed */
      if (hc_final * hc_init < 0.0) {
        T_test = zbrent(heatcool, fmin(T_f, T), fmax(T_f, T), 1.0);
        T_f = T_test;
      }
    }

    /* Apply temperature floor */
    T_f = MAX(T_f, g_minCoolingTemp);

    /* Convert final temperature back to internal energy */
    E_f = T_f / (2.0 / 3.0) * (n * CONST_kB);

    /* Convert back to pressure in code units */
    p_f = E_f * (g_gamma - 1) / UNIT_PRESSURE;

    /* Update the pressure in the cell */
    VV[PRS][k][j][i] = p_f;
  }
}

/*!
 * \brief Computes the net heating or cooling rate at temperature T.
 *
 * \param [in] T   Temperature in Kelvin.
 * \return     Net heating/cooling rate.
 */
double heatcool(double T)
{

  double lambda;                      // Net heating/cooling rate
  double sqT = sqrt(T);
  double comp_heat, comp_cool, xray_heat;
  double line_cool, brem_cool;

  #if LINE_DRIVEN_WIND != NO
    /* Compute electron number density using ne_rat function */
    ne = nH * ne_rat(T);
  #else
    /* Assume fully ionized gas */
    ne = 1.21 * nH;
  #endif

  /* Heating and cooling rates with precomputed coefficients */
  comp_heat = comp_h_pre * (8.9e-36 * xi * tx);
  comp_cool = comp_c_pre * (8.9e-36 * xi * (4.0 * T));
  xray_heat = xray_h_pre * (1.5e-21 * (sqsqxi / sqT));
  line_cool = line_c_pre * ((1e-16 * exp(-1.3e5 / T) / sqxi / T) +
           fmin(fmin(1e-24, 5e-27 * sqT), 1.5e-17 / T));
  brem_cool = brem_c_pre * (3.3e-27 * sqT);

  /* Net heating/cooling rate */
  lambda = nH * (ne * comp_heat + nH * xray_heat
              -  ne * comp_cool - ne * line_cool - ne * brem_cool);

  return lambda;
}

/*!
 * \brief Computes the electron number density ratio ne/nH at given temperature.
 *
 * \param [in] temp  Temperature in Kelvin.
 * \return     Electron number density ratio ne/nH.
 */
double ne_rat(double T)
{
  double ne_ratio;
  if (T < 1.5e4){
    ne_ratio = 1e-2 + pow(10, (-51.59417133 + 12.27740153 * log10(T)));
  }else if (T >= 1.5e4 && T < 3.3e4){
    ne_ratio = pow(10, (-3.80749689 + 0.86092628 * log10(T)));
  } else{
    ne_ratio = 1.21;
  }

  return ne_ratio;
}

/*!
 * \brief Function used in root finding to solve for temperature.
 *
 * Calculates the difference between the internal energy at temperature `temp`,
 * minus the current energy `E`, minus the average heating/cooling rate over the time step `dt_share`.
 *
 * \param [in] temp  Temperature in Kelvin.
 * \return     Value of the function at `temp`.
 */
double zfunc(double temp)
{
  double ans = (temp * n * CONST_kB / (2.0 / 3.0)) - E - dt_share * (hc_init + heatcool(temp)) / 2.0;
  return ans;
}

/*!
 * \brief Root-finding function using Brent's method.
 *
 * Finds the root of a function `func` in the interval [x1, x2] to within a specified tolerance `tol`.
 *
 * \param [in] func  Function pointer to the function whose root is sought.
 * \param [in] x1    Lower bound of the interval.
 * \param [in] x2    Upper bound of the interval.
 * \param [in] tol   Desired tolerance.
 * \return     Root of the function `func` within the interval [x1, x2].
 */
double zbrent(double (*func)(double), double x1, double x2, double tol)
{
  int iter;
  double a = x1, b = x2, c = x2;  // Initialization
  double d = 0.0, e = 0.0;        // Interval variables
  double fa = func(a), fb = func(b), fc = fb;  // Function evaluations
  double p, q, r, s, tol1, xm;    // Algorithm variables

  if (fb * fa > 0.0) {
    printf("ZBRENT: Function must bracket zero. f(%e)=%e, f(%e)=%e\n", x1, fa, x2, fb);
    return b;
  }

  for (iter = 1; iter <= ITMAX; iter++) {
    if (fb * fc > 0.0) {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
    xm = 0.5 * (c - b);

    if (fabs(xm) <= tol1 || fb == 0.0)
      return b;

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = fabs(p);
      double min1 = 3.0 * xm * q - fabs(tol1 * q);
      double min2 = fabs(e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    fb = func(b);
  }
  printf("Maximum number of iterations exceeded in zbrent\n");
  return b;
}
