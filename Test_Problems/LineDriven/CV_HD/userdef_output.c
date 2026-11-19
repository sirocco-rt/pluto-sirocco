#include "pluto.h"

/* ====================================================================
   Accessor functions defined at the bottom of line_connect.c.
   Declared here with 'extern' so no header changes are needed.
   ==================================================================== */
extern int    LC_GetNFluxAngles (void);
extern int    LC_GetMPoints     (void);
extern int    LC_IsMUVAllZero   (int k, int j, int i);
extern double LC_GetFluxR       (int a, int k, int j, int i);
extern double LC_GetFluxT       (int a, int k, int j, int i);
extern double LC_GetDvDs        (int a, int k, int j, int i);
extern double LC_GetMUVFit      (int m, int k, int j, int i);
extern double LC_GetTFit        (int m);

/* *************************************************************** */
/*!
 * \brief Defines user-defined output variables for the PLUTO simulation.
 *
 * Computes temperature, ionization parameter, number densities, force
 * multipliers M_max1 / M_max2 (legacy M_UV variable name; cak_v4 includes all eligible
 * line transitions across the supplied spectral range),
 * the Sobolev optical depth t, and the velocity gradient |dv/ds|.
 *
 * M_max1  : maximum M_UV across all angular flux bins
 * M_max2  : multiplier at the angle of maximum SIROCCO UV-band flux magnitude
 *           (the directional-flux interface remains UV-band only)
 * t       : Sobolev optical depth t = sigma_e * rho * v_th / |dv/ds| (all-frequency CAK multiplier table)
 * dv_ds   : velocity gradient |dv/ds| at the dominant angle
 *
 * \authors Nick Higginbottom, Amin Mosallanezhad
 */
void ComputeUserVar(const Data *d, Grid *grid)
{
  int i, j, k;
  int nv;
  int ii;
  int iangle;

  #if COOLING != NO
    double ***comp_heat, ***comp_cool;
    double ***xray_heat, ***line_cool, ***brem_cool;
    double ***comp_heat_pre, ***comp_cool_pre;
    double ***line_cool_pre, ***brem_cool_pre;
    double ***xray_heat_pre, ***xi_ion_pre;
  #endif

  double ***T_out, ***Tr_out, ***xi_out;
  double ***ne_out, ***nH_out;
  double t_UV, M_UV, t_UV2, dvds2;

  double ***gr_out, ***gt_out, ***gp_out;
  double grad[3];
  double Vc[NVAR];

  double ***dv_ds_out;
  double ***t_out;
  double ***M_out1, ***M_out2;
  double M_max, eta_max, M_temp1, M_temp2, fmax;

  double rho, p, T, xi, ne, nH, n;
  double mu, lx, tx, r;
  double sigma_e, v_th;
  double sqrt_T, sqxi, sqsqxi;

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  int    NFLUX_ANGLES, MPOINTS;
  double flux_r_val, flux_t_val, dvds_val;
  double *M_UV_array = NULL;

  #if COOLING != NO
    comp_heat      = GetUserVar("ch");
    comp_cool      = GetUserVar("cc");
    xray_heat      = GetUserVar("xh");
    line_cool      = GetUserVar("lc");
    brem_cool      = GetUserVar("bc");
    comp_heat_pre  = GetUserVar("ch_pre");
    comp_cool_pre  = GetUserVar("cc_pre");
    xray_heat_pre  = GetUserVar("xh_pre");
    line_cool_pre  = GetUserVar("lc_pre");
    brem_cool_pre  = GetUserVar("bc_pre");
    xi_ion_pre     = GetUserVar("XI_pre");
  #endif

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

  mu      = g_inputParam[MU];
  lx      = g_inputParam[L_STAR] * g_inputParam[F_X];
  /* Must match LineForce() and cak_v4 (fixed sigma_T/m_p). */
  sigma_e = 0.66524e-24 / 1.672661e-24;

  /* Fetch lc_ dimensions once per timestep */
  NFLUX_ANGLES = LC_GetNFluxAngles();
  MPOINTS      = LC_GetMPoints();

  /* Allocate the per-cell all-frequency CAK M(t) lookup table.  The table
   retains its historical M_UV_data.dat / M_UV symbol names. */
  if (MPOINTS > 0) {
    M_UV_array = (double *) malloc(MPOINTS * sizeof(double));
  }

  DOM_LOOP(k, j, i) {

    r   = grid->x[IDIR][i] * UNIT_LENGTH;
    rho = d->Vc[RHO][k][j][i] * UNIT_DENSITY;

    /* Temperature */
    #if EOS == ISOTHERMAL
      T = T_out[k][j][i] = g_inputParam[T_ISO];
    #else
      T = T_out[k][j][i] = d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i] * KELVIN * mu;
    #endif

    /* Hydrogen number density */
    nH = rho / (1.43 * CONST_mp);
    nH_out[k][j][i] = nH;

    n = rho / (mu * CONST_mp);

    /* Electron number density */
    double ne_ratio;
    if (T < 1.5e4) {
      ne_ratio = 1e-2 + pow(10.0, (-51.59417133 + 12.27740153 * log10(T)));
    } else if (T < 3.3e4) {
      ne_ratio = pow(10.0, (-3.80749689 + 0.86092628 * log10(T)));
    } else {
      ne_ratio = 1.21;
    }

    #if COOLING != NO
      ne = nH * ne_ratio;
    #else
      ne = 1.21 * nH;
    #endif

    ne_out[k][j][i] = ne;

    /* Thermal velocity (cm/s) */
    v_th = sqrt(2.0 * CONST_kB * T / CONST_mp);

    #if COOLING != NO
      /* Ionization parameter and radiation temperature.
       * Use analytical estimate before flux files are loaded (NFLUX_ANGLES == 0,
       * i.e. before cycle 0 completes).  Use SIROCCO output once available.
       * No sentinel check needed: read_sirocco_heatcool flag=0 now writes
       * real analytical xi, so sirocco_xi is never negative. */
      if (NFLUX_ANGLES == 0) {
        xi = lx / (nH * r * r);
        xi_out[k][j][i] = xi;
        Tr_out[k][j][i] = g_inputParam[T_X];
      } else {
        xi = d->sirocco_xi[k][j][i];
        xi_out[k][j][i] = xi;
        Tr_out[k][j][i] = d->sirocco_t_r[k][j][i];
      }
      tx = Tr_out[k][j][i];

      sqrt_T = sqrt(T);
      sqxi   = sqrt(xi);
      sqsqxi = pow(xi, 0.25);

      /* Precomputed prefactors */
      comp_heat_pre[k][j][i] = d->comp_h_pre[k][j][i];
      comp_cool_pre[k][j][i] = d->comp_c_pre[k][j][i];
      xray_heat_pre[k][j][i] = d->xray_h_pre[k][j][i];
      line_cool_pre[k][j][i] = d->line_c_pre[k][j][i];
      brem_cool_pre[k][j][i] = d->brem_c_pre[k][j][i];
      xi_ion_pre[k][j][i]    = d->xi_ion_pre[k][j][i];

      /* Apply prefactors only once SIROCCO output is available.
       * Before flux files load (NFLUX_ANGLES == 0), use bare analytical rates. */
      if (NFLUX_ANGLES == 0) {
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
    #endif  /* COOLING != NO */


    #if LINE_DRIVEN_WIND != NO

      /* Body force acceleration via LineForce() */
      NVAR_LOOP(nv) Vc[nv] = d->Vc[nv][k][j][i];
      LineForce(Vc, grad, x1[i], x2[j], x3[k], i, j, k);

      #if (BODY_FORCE & VECTOR)
        gr_out[k][j][i] = grad[0];
        gt_out[k][j][i] = grad[1];
        gp_out[k][j][i] = grad[2];
      #endif

      /* Maximum allowed force multiplier (CAK theory) */
      M_max   = 4400.0;
      eta_max = 7.954346e+07;

      /* Recompute temperature for EOS consistency */
      #if EOS != ISOTHERMAL
        T = Vc[PRS] / Vc[RHO] * KELVIN * mu;
      #else
        T = g_inputParam[T_ISO];
      #endif

      /* Recompute thermal velocity with current cell T */
      v_th = sqrt(2.0 * CONST_kB * T / CONST_mp);

      /* Populate per-cell M_UV lookup array from CAK table */
      if (MPOINTS > 0 && M_UV_array != NULL) {
        for (ii = 0; ii < MPOINTS; ii++) {
          M_UV_array[ii] = LC_GetMUVFit(ii, k, j, i);
        }
      }

      /* -----------------------------------------------------------
         Per-angle loop: compute M_UV (Sobolev force multiplier),
         t_UV (Sobolev optical depth), and dv/ds for each angular bin.

         M_temp1 = maximum M_UV across all angles            -> M_max1
         M_temp2 = M_UV at angle of maximum flux magnitude   -> M_max2
         t_UV2   = Sobolev optical depth at dominant angle   -> t
         dvds2   = |dv/ds| at dominant angle                 -> dv_ds
         ----------------------------------------------------------- */
      M_temp1 = -1.0;
      M_temp2 = -1.0;
      fmax    = 1e-99;
      t_UV2   = 0.0;
      dvds2   = 0.0;

      for (iangle = 0; iangle < NFLUX_ANGLES; iangle++) {

        dvds_val   = LC_GetDvDs (iangle, k, j, i);
        flux_r_val = LC_GetFluxR(iangle, k, j, i);
        flux_t_val = LC_GetFluxT(iangle, k, j, i);

        if (dvds_val > 0.0) {
          /* Sobolev optical depth for this angle */
          t_UV = sigma_e * rho * v_th / dvds_val;

          /* Force multiplier M_UV from CAK table or zero */
          if (MPOINTS > 0 && M_UV_array != NULL &&
              !LC_IsMUVAllZero(k, j, i)) {
            /* Linear interpolation in log10(t) – log10(M) space */
            double log_t = log10(t_UV);
            double log_M;
            double t0v = LC_GetTFit(0);
            double t1v = LC_GetTFit(MPOINTS - 1);
            if (log_t <= t0v) {
              log_M = M_UV_array[0];
            } else if (log_t >= t1v) {
              log_M = M_UV_array[MPOINTS - 1];
            } else {
              int m;
              log_M = M_UV_array[MPOINTS - 1];  /* fallback */
              for (m = 1; m < MPOINTS; m++) {
                t1v = LC_GetTFit(m);
                t0v = LC_GetTFit(m - 1);
                if (log_t < t1v) {
                  double frac = (log_t - t0v) / (t1v - t0v);
                  log_M = M_UV_array[m-1] + frac * (M_UV_array[m] - M_UV_array[m-1]);
                  break;
                }
              }
            }
            M_UV = pow(10.0, log_M);
          } else {
            M_UV = 0.0;   /* no CAK table loaded yet */
          }

          if (M_UV > M_max) M_UV = M_max;

          /* Track maximum M_UV across all angles */
          if (M_UV > M_temp1) M_temp1 = M_UV;

        } else {
          M_UV = 0.0;
          t_UV = 0.0;
        }

        /* Track angle with maximum flux magnitude for M_temp2 */
        double flux_sq = flux_r_val * flux_r_val + flux_t_val * flux_t_val;
        if (flux_sq > fmax) {
          fmax    = flux_sq;
          M_temp2 = M_UV;
          t_UV2   = t_UV;
          dvds2   = dvds_val;
        }

      }  /* end iangle loop */

      /* Output:
         M_max1 = maximum M_UV (force multiplier) across all angles
         M_max2 = M_UV at angle of maximum UV flux magnitude
         t      = Sobolev optical depth at dominant angle
         dv_ds  = |dv/ds| at dominant angle                         */
      M_out1[k][j][i]    = (M_temp1 >= 0.0) ? M_temp1 : 0.0;
      M_out2[k][j][i]    = (M_temp2 >= 0.0) ? M_temp2 : 0.0;
      t_out[k][j][i]     = t_UV2;
      dv_ds_out[k][j][i] = dvds2;

    #endif  /* LINE_DRIVEN_WIND != NO */

  }  /* End DOM_LOOP */

  if (M_UV_array) { free(M_UV_array); M_UV_array = NULL; }
}

/* **************************************************************** */
void ChangeOutputVar()
{
  Image *image;

  #ifdef PARTICLES
  /* SetOutputVar("energy", PARTICLES_FLT_OUTPUT, NO); */
  #endif
}
