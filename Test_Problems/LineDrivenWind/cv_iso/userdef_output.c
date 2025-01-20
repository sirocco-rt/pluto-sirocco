#include "pluto.h"

/* *************************************************************** */
/*!
 * \brief Defines user-defined output variables for the PLUTO simulation.
 *
 * This function computes various physical quantities such as temperature,
 * ionization parameter, number densities, and force multipliers. These
 * quantities are stored in user-defined variables for output and analysis.
 *
 * \param [in]  d     Pointer to the Data structure containing simulation variables.
 * \param [in]  grid  Pointer to the Grid structure.

 \authors Nick Higginbottom
 \date    July 28, 2018
 revised by: Amin Mosallanezhad

 */
void ComputeUserVar(const Data *d, Grid *grid)

{

  int i, j, k;                   // Loop indices for grid traversal
  int nv;                        // Variable index for NVAR_LOOP
  int ii;                        // General-purpose index
  int iangle;                    // Index for angular bins in radiation flux

  #if COOLING != NO
    double ***comp_heat, ***comp_cool;
    double ***xray_heat, ***line_cool, ***brem_cool;
    //
    double ***comp_heat_pre, ***comp_cool_pre;
    double ***line_cool_pre, ***brem_cool_pre;
    double ***xray_heat_pre, ***xi_ion_pre;
  #endif

  double ***T_out, ***Tr_out, ***xi_out;     // Temperature and ionization parameter
  double ***ne_out, ***nH_out;               // Electron and hydrogen number densities
  double krad, alpharad;          // Radiation force parameters
  double t_UV, M_UV, t_UV2, dvds2;              // Optical depth parameter and force multiplier

  #if LINE_DRIVEN_WIND != NO
    double M_UV_array[MPOINTS];     // Array for force multipliers from fits
  #endif

  // #if (BODY_FORCE & VECTOR)
    double ***gr_out, ***gt_out, ***gp_out;
    double grad[3];                 // Acceleration vector [grad_r, grad_t, grad_p]
    double Vc[NVAR];                // Array for cell variables
  // #endif

  double ***dv_ds_out;            // Output array for velocity gradients
  double ***t_out;                // Output array for optical depth parameter
  double ***M_out1, ***M_out2;    // Output arrays for force multipliers
  double M_max, M_min, eta_max, M_temp1, M_temp2, fmax;

  // Physical variables
  double rho, p, T, xi, ne, nH, n;  // Density, pressure, temperature, ionization parameter, number densities
  double mu, lx, tx, r;             // Mean molecular weight, luminosity, temperature, radius
  double sigma_e, v_th;             // Electron scattering opacity, thermal velocity
  double sqrt_T, sqxi, sqsqxi;
//
  double *x1 = grid->x[IDIR];       // Radial coordinate array
  double *x2 = grid->x[JDIR];       // Angular coordinate array
  double *x3 = grid->x[KDIR];       // Azimuthal coordinate array (if applicable)

  #if COOLING != NO
    // Get pointers to user-defined variables for heating and cooling rates
    comp_heat = GetUserVar("ch");
    comp_cool = GetUserVar("cc");
    xray_heat = GetUserVar("xh");
    line_cool = GetUserVar("lc");
    brem_cool = GetUserVar("bc");
//
    comp_heat_pre  = GetUserVar("ch_pre");
    comp_cool_pre  = GetUserVar("cc_pre");
    xray_heat_pre  = GetUserVar("xh_pre");
    line_cool_pre  = GetUserVar("lc_pre");
    brem_cool_pre  = GetUserVar("bc_pre");
    xi_ion_pre     = GetUserVar("XI_pre");
  #endif

  // Get pointers to user-defined variables for output
  ne_out  = GetUserVar("ne");
  nH_out  = GetUserVar("nh");
  T_out   = GetUserVar("T");
  xi_out  = GetUserVar("XI");
  Tr_out  = GetUserVar("T_r");

  #if (BODY_FORCE & VECTOR)
    gr_out  = GetUserVar("gr");
    gt_out  = GetUserVar("gt");
    gp_out  = GetUserVar("gp");
  #endif

  dv_ds_out = GetUserVar("dv_ds");
  t_out     = GetUserVar("t");
  M_out1    = GetUserVar("M_max1");
  M_out2    = GetUserVar("M_max2");

  // Retrieve input parameters from global input parameter array
  mu = g_inputParam[MU];                               // Mean molecular weight
  lx = g_inputParam[L_star] * g_inputParam[f_x];       // X-ray luminosity (erg/s)
  tx = g_inputParam[T_x];                              // X-ray temperature (K)

  sigma_e = CONST_sigmaT / (CONST_amu * 1.18);         // Effective electron scattering opacity (cm^2/g)


  // Loop over the entire computational domain
  DOM_LOOP(k, j, i) {
    // Compute radius in physical units
    r = grid->x[IDIR][i] * UNIT_LENGTH;

    // Compute density in physical units (g/cm^3)
    rho = d->Vc[RHO][k][j][i] * UNIT_DENSITY;

    // Compute temperature (K)
    #if EOS == ISOTHERMAL
      T = T_out[k][j][i] = g_inputParam[T_ISO];          // Isothermal temperature
    #else
      T = T_out[k][j][i] = d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i] * KELVIN * mu;
    #endif

    // Compute hydrogen number density (cm^-3) assuming stellar abundances
    nH = rho / (1.43 * CONST_mp);
    nH_out[k][j][i] = nH;


    // Compute particle number density (cm^-3)
    n = rho / (mu * CONST_mp);

    // Electron number density assuming full ionization (cm^-3)
    double ne_ratio;
    if (T < 1.5e4){
      ne_ratio = 1e-2 + pow(10, (-51.59417133 + 12.27740153 * log10(T)));
    }else if (T >= 1.5e4 && T < 3.3e4){
      ne_ratio = pow(10, (-3.80749689 + 0.86092628 * log10(T)));
    } else{
      ne_ratio = 1.21;
    }

    #if COOLING != NO
      ne = nH * ne_ratio;
    #else
      ne = 1.21 * nH;
    #endif

    ne_out[k][j][i] = ne;

    // Compute thermal velocity for hydrogen (cm/s)
    v_th = sqrt(2.0 * CONST_kB * T / CONST_mp);

    #if COOLING != NO
      // Compute ionization parameter xi (erg cm/s) and radiation Temperature
      if (g_time <= 1.0) {
        xi = lx / (nH * r * r);             // Initial ionization parameter
        xi_out[k][j][i] = xi;
        Tr_out[k][j][i] = g_inputParam[T_x];
      } else {
        // xi = sirocco_xi[k][j][i] * xi_ion_pre[k][j][i];            // Updated ionization parameter from SIROCCO
        xi = d->sirocco_xi[k][j][i];
        xi_out[k][j][i] = xi;
        Tr_out[k][j][i] = d->sirocco_t_r[k][j][i];
      }
//
      sqrt_T = sqrt(T);
      sqxi   = sqrt(xi);           // xi^(0.5)
      sqsqxi = pow(xi, 0.25);      // xi^(0.25)
//
      // Store precomputed heating and cooling rates
      comp_heat_pre[k][j][i] = d->comp_h_pre[k][j][i];
      comp_cool_pre[k][j][i] = d->comp_c_pre[k][j][i];
      xray_heat_pre[k][j][i] = d->xray_h_pre[k][j][i];
      line_cool_pre[k][j][i] = d->line_c_pre[k][j][i];
      brem_cool_pre[k][j][i] = d->brem_c_pre[k][j][i];
      xi_ion_pre[k][j][i]    = d->xi_ion_pre[k][j][i];

      if (g_time < 3.0) {
        comp_heat[k][j][i] = (8.9e-36 * xi * tx);
        comp_cool[k][j][i] = (8.9e-36 * xi * (4.0 * T));
        xray_heat[k][j][i] = (1.5e-21 * (sqsqxi / sqrt_T));
        line_cool[k][j][i] = ((1e-16 * exp(-1.3e5 / T) / sqxi / T) +
                             fmin(fmin(1e-24, 5e-27 * sqrt_T), 1.5e-17 / T));
        brem_cool[k][j][i] = (3.3e-27 * sqrt_T);
      } else {
        comp_heat[k][j][i] = comp_heat_pre[k][j][i] * (8.9e-36 * xi * tx);
        comp_cool[k][j][i] = comp_cool_pre[k][j][i] * (8.9e-36 * xi * (4.0 * T));
        xray_heat[k][j][i] = xray_heat_pre[k][j][i] * (1.5e-21 * (sqsqxi / sqrt_T));
        line_cool[k][j][i] = line_cool_pre[k][j][i] * ((1e-16 * exp(-1.3e5 / T) / sqxi / T) +
                             fmin(fmin(1e-24, 5e-27 * sqrt_T), 1.5e-17 / T));
        brem_cool[k][j][i] = brem_cool_pre[k][j][i] * (3.3e-27 * sqrt_T);
      }

    #endif


      #if LINE_DRIVEN_WIND != NO

      // Copy cell variables into Vc array for body force calculation
      NVAR_LOOP(nv)  Vc[nv] = d->Vc[nv][k][j][i];


      // Compute body force vector at the current cell position

      LineForce(Vc, grad, x1[i], x2[j], x3[k], i, j, k);
      //
      // // Store acceleration components in spherical coordinates
      gr_out[k][j][i] = grad[0];
      gt_out[k][j][i] = grad[1];
      gp_out[k][j][i] = grad[2];

      // Maximum allowed force multiplier (from CAK theory)
      M_max = 4400.;
      eta_max = 7.954346e+07;

      krad = g_inputParam[KRAD];
      alpharad = g_inputParam[ALPHARAD];

      // Compute temperature for both isothermal and ideal equation of state
      #if EOS != ISOTHERMAL
        T = Vc[PRS] / Vc[RHO] * KELVIN * mu;
      #else
        T = g_inputParam[T_ISO];
      #endif

    // If krad and alpharad are special values (999), use precomputed force multipliers
    if (krad == 999 && alpharad == 999) {
      // Copy M_UV_fit values into M_UV_array for the current cell
      for (ii = 0; ii < MPOINTS; ii++) {
        M_UV_array[ii] = M_UV_fit[ii][k][j][i];
      }
    }

    M_temp1 = -1.0;  // Temporary variable to store maximum M_UV
    M_temp2 = -1.0;  // Temporary variable for flux-weighted M_UV
    fmax    = 1e-99; // Maximum flux magnitude (initially very small)

    for (iangle = 0; iangle < NFLUX_ANGLES; iangle++) {
//
  	    if (dvds_array[iangle][k][j][i] > 0.0) {
  	        t_UV = sigma_e * rho * v_th / dvds_array[iangle][k][j][i];
  		      if (krad == 999 && alpharad == 999) {
  				     M_UV = linterp(log10(t_UV), t_fit, M_UV_array, MPOINTS);
  			    } else {
  		         M_UV = krad * pow(t_UV, alpharad);
  			    }

  	        if (M_UV > M_max)   M_UV = M_max;
  			    if (M_UV > M_temp1) M_temp1 = M_UV;
  	      }

       if (pow(flux_r_UV[iangle][k][j][i], 2) + pow(flux_t_UV[iangle][k][j][i], 2) > fmax) {
      		fmax = pow(M_UV * flux_r_UV[iangle][k][j][i], 2) + pow(M_UV * flux_t_UV[iangle][k][j][i], 2);
      		M_temp2 = M_UV;

          t_UV2   = t_UV;
          dvds2   = dvds_array[iangle][k][j][i];
       }
//
    }

    // Store maximum force multipliers in output arrays
    M_out1[k][j][i] = M_temp1;  // Maximum M_UV across all angles
    M_out2[k][j][i] = M_temp2;  // M_UV corresponding to maximum flux

    dv_ds_out[k][j][i] = dvds2;
    t_out[k][j][i]     = t_UV2;


    #endif


  }  // End of DOM_LOOP

}

/* ************************************************************* */
/*!
 * \brief Allows modification of output variables.
 *
 * This function can be customized to enable or disable the output of specific variables
 * in the simulation data. Currently, it is empty but provides a template for future use.
 */
void ChangeOutputVar()
{
  Image *image;  // Placeholder for image handling (unused here)

  #ifdef PARTICLES
  // Example: Uncomment and modify the following lines to change particle output variables
  // SetOutputVar("energy", PARTICLES_FLT_OUTPUT, NO);
  // SetOutputVar("x1",     PARTICLES_FLT_OUTPUT, NO);
  // SetOutputVar("vx1",    PARTICLES_FLT_OUTPUT, NO);
  #endif
}
