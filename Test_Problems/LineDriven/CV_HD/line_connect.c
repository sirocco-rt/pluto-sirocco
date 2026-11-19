/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  PLUTO-SIROCCO line-driven wind connection module (cv_ideal).

  Self-contained module: all state is declared as file-scope static
  variables.

  Public API (declared in Src/prototypes.h):
    void read_sirocco_fluxes   (Data *, Grid *)
    void read_sirocco_heatcool (Data *, Grid *, int)
    void VGradCalc             (const Data *, Grid *)
    void LineForce             (double *, double *, double, double,
                                double, int, int, int)

  Private helpers (static, not in any header):
    static double lc_linterp  (double, double *, double *, int)
    static void   lc_bilinear (double[2], double[2], double[2],
                                double[2], double[2], double[2],
                                double[2], double[2], int)

  \author  A. Mosallanezhad  (a.mosallanezhad@soton.ac.uk)
  \date    April 2026
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* ====================================================================
   Module-level state — all static so nothing leaks into the global
   namespace and no extern declarations are needed in pluto.h /
   globals.h.
   ==================================================================== */

#define LC_LINELENGTH 800  /* max chars per input-file line */

/* ---- scalar counters ---- */
static int lc_NFLUX_ANGLES = 0;  /*!< Number of angular flux bins (e.g. 36)  */
static int lc_MPOINTS      = 0;  /*!< Number of force-multiplier table rows  */

/* ---- 4-D flux arrays [angle][k][j][i] ---- */
static double ****lc_flux_r = NULL; /*!< Radial UV flux per angular bin       */
static double ****lc_flux_t = NULL; /*!< Theta  UV flux per angular bin       */
static double ****lc_flux_p = NULL; /*!< Phi    UV flux per angular bin       */

/* ---- force-multiplier table ---- */
static double ****lc_M_UV_fit = NULL; /*!< log10(M) per cell [mpt][k][j][i]  */
static double  *lc_t_fit      = NULL; /*!< log10(t) knot values [mpt]         */
static int    ***lc_M_UV_all_zero = NULL; /*!< 1 when every M(t) knot is <= 0 */

/* ---- velocity-gradient workspace ---- */
static double ****lc_dvds_r_off   = NULL; /*!< r-offset for dv/ds step        */
static double ****lc_dvds_t_off   = NULL; /*!< theta-offset for dv/ds step    */
static double ****lc_dvds_mod_off = NULL; /*!< path length ds for dv/ds step  */
static double ****lc_dvds_array   = NULL; /*!< |dv/ds| per angle per cell     */

/* ---- interpolation corner locations [2][k][j][i] ---- */
static double ****lc_x11 = NULL; /*!< lower-left vertex (r, theta) in CGS    */
static double ****lc_x22 = NULL; /*!< upper-right vertex (r, theta) in CGS   */

/* ---- first-call guard for VGradCalc geometry setup ---- */
static int lc_dvds_ready = 0; /*!< 0 = geometry not yet built; 1 = done      */

/* ---- private helper forward declarations ---- */
static double lc_linterp  (double x, double *xarr, double *yarr, int n);
static void   lc_bilinear (double x11[2], double x22[2],
                            double v11[2], double v12[2],
                            double v21[2], double v22[2],
                            double test[2], double ans[2], int flag);


/* ====================================================================
   1.  read_sirocco_fluxes
       Reads directional_flux_{r,theta,phi}.dat and the all-frequency M_UV_data.dat from cak_v4.
   ==================================================================== */
void read_sirocco_fluxes (Data *d, Grid *grid)
{
    int    i, j, iflux, iaxis, icount, match, m;
    long   int k, ii, jj, kk, nwords;
    FILE   *fptr;
    char   aline[LC_LINELENGTH];
    double x1in, x2in, temp;
    double tol = 1e-6;
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];

    print ("  > read_sirocco_fluxes: IBEG=%d JBEG=%d\n", IBEG, JBEG);

    /* ---- read flux files for r, theta, phi ---- */
    for (iaxis = 0; iaxis < 3; iaxis++) {

        const char *fname = (iaxis == 0) ? "directional_flux_r.dat"
                          : (iaxis == 1) ? "directional_flux_theta.dat"
                                         : "directional_flux_phi.dat";

        fptr = fopen (fname, "r");
        if (fptr == NULL) {
            print ("  > read_sirocco_fluxes: no file %s\n", fname);
            continue;
        }

        /* skip first header line, read second (contains NANGLES) */
        if (fgets (aline, LC_LINELENGTH, fptr) == NULL ||
            fgets (aline, LC_LINELENGTH, fptr) == NULL) {
            print ("! read_sirocco_fluxes: header read error in %s\n", fname);
            QUIT_PLUTO(1);
        }

        nwords = sscanf (aline, "%*s %*s %ld", &ii);
        if (nwords != 1) {
            print ("! read_sirocco_fluxes: header format error in %s\n", fname);
            QUIT_PLUTO(1);
        }

        if (iaxis == 0) {
            lc_NFLUX_ANGLES = (int)ii;
            print ("  > %d angular flux bins\n", lc_NFLUX_ANGLES);
            lc_flux_r = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);
        } else if (iaxis == 1) {
            if ((int)ii != lc_NFLUX_ANGLES) {
                print ("! read_sirocco_fluxes: NANGLES mismatch in theta file\n");
                QUIT_PLUTO(1);
            }
            lc_flux_t = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);
        } else {
            if ((int)ii != lc_NFLUX_ANGLES) {
                print ("! read_sirocco_fluxes: NANGLES mismatch in phi file\n");
                QUIT_PLUTO(1);
            }
            lc_flux_p = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);
        }

        icount = match = 0;

        while (fscanf (fptr, "%ld ", &ii) != EOF) {
            if (fscanf (fptr, "%ld %ld %le %le", &jj, &kk, &x1in, &x2in) != 4) {
                print ("! read_sirocco_fluxes: data read error in %s\n", fname);
                QUIT_PLUTO(1);
            }

            DOM_LOOP (k, j, i) {
                if (fabs (1.0 - x1in / UNIT_LENGTH / x1[i]) < tol &&
                    fabs (1.0 - x2in / x2[j])              < tol) {
                    for (iflux = 0; iflux < lc_NFLUX_ANGLES; iflux++) {
                        if (fscanf (fptr, "%le", &temp) != 1) {
                            print ("! read_sirocco_fluxes: flux read error\n");
                            QUIT_PLUTO(1);
                        }
                        if (iaxis == 0) lc_flux_r[iflux][k][j][i] = temp;
                        if (iaxis == 1) lc_flux_t[iflux][k][j][i] = temp;
                        if (iaxis == 2) lc_flux_p[iflux][k][j][i] = temp;
                    }
                    match = 1;
                    icount++;
                }
            }

            if (!match) {
                for (iflux = 0; iflux < lc_NFLUX_ANGLES; iflux++) {
                    if (fscanf (fptr, "%le", &temp) != 1) {
                        print ("! read_sirocco_fluxes: skip-read error\n");
                        QUIT_PLUTO(1);
                    }
                }
            }
            match = 0;
        }

        {
            /* Sum icount across all MPI processes to report the global
             * total alongside the per-process count.                    */
            int icount_global = icount;
#ifdef PARALLEL
            MPI_Allreduce (&icount, &icount_global, 1, MPI_INT, MPI_SUM,
                           MPI_COMM_WORLD);
#endif
            print ("  > read %d cells from %s (this proc; %d total)\n",
                   icount, fname, icount_global);
        }
        fclose (fptr);
    }

    /* ---- force-multiplier table (loaded whenever present) ---- */
    {

        fptr = fopen ("M_UV_data.dat", "r");
        if (fptr == NULL) {
            print ("  > read_sirocco_fluxes: no M_UV_data.dat yet, using pure electron scattering\n");
            /* lc_MPOINTS stays 0 --> LineForce uses M_UV=0 */
        } else {

        if (fgets (aline, LC_LINELENGTH, fptr) == NULL) {
            print ("! read_sirocco_fluxes: M_UV header read error\n");
            QUIT_PLUTO(1);
        }

        nwords = sscanf (aline, "%*s %ld", &ii);
        if (nwords != 1) {
            print ("! read_sirocco_fluxes: M_UV header format error\n");
            QUIT_PLUTO(1);
        }
        lc_MPOINTS = (int)ii;
        print ("  > %d force-multiplier table rows\n", lc_MPOINTS);

        /* ---- Read t_fit knots, skipping any bad (non-positive) entries ----
         *
         * Clamping bad t values to a tiny positive number (e.g. 1e-20) and the
         * corresponding M values to 1e-30 creates a spurious anchor at the far
         * left of the table with log10(M) = -30.  Because M must decrease as t
         * increases (M = k*t^{-alpha}), this anchor has the wrong sign of slope
         * relative to all subsequent valid knots, inverting the interpolation in
         * the optically-thin regime and returning near-zero M where it should be
         * large.  The monotonicity check does not catch this because -20 < -3 etc.
         *
         * The correct fix is to skip bad rows entirely so the table always has a
         * physically meaningful, monotone-increasing set of log10(t) knots with a
         * consistent (negative) slope in log10(M).  A temporary raw-value array
         * is used to collect valid entries; lc_MPOINTS is then updated to the
         * number of valid knots before the M_UV_fit 4-D array is allocated.     */

        {
            int    n_raw = lc_MPOINTS;   /* total rows declared in file header */
            int    n_valid = 0;          /* rows that survive the validity test */
            double *t_raw;               /* temporary storage for raw t values  */
            int   *valid_idx;            /* which file-row indices are valid     */

            t_raw     = (double *) calloc (n_raw, sizeof (double));
            valid_idx = (int    *) calloc (n_raw, sizeof (int));

            /* skip one token (header label for t values) */
            if (fscanf (fptr, "%*s") == EOF) { /* intentionally empty */ }

            for (iflux = 0; iflux < n_raw; iflux++) {
                if (fscanf (fptr, "%le", &temp) != 1) {
                    print ("! read_sirocco_fluxes: t_fit read error at index %d\n", iflux);
                    QUIT_PLUTO(1);
                }
                if (temp > 0.0 && isfinite (temp)) {
                    t_raw[n_valid]     = temp;
                    valid_idx[n_valid] = iflux;  /* original column index in file */
                    n_valid++;
                } else {
                    print ("  > WARNING: skipping non-positive t_fit[%d]=%e "
                           "(would invert M(t) slope if clamped)\n", iflux, temp);
                }
            }

            /* Update lc_MPOINTS to the number of usable knots */
            lc_MPOINTS = n_valid;
            print ("  > %d valid t_fit knots (%d total, %d skipped)\n",
                   n_valid, n_raw, n_raw - n_valid);

            if (lc_MPOINTS == 0) {
                print ("  > WARNING: no valid t_fit knots — disabling CAK table\n");
                free (t_raw); free (valid_idx);
                lc_t_fit = NULL;
                fclose (fptr);
                goto end_muv_block;
            }

            /* Build lc_t_fit from valid entries */
            lc_t_fit = (double *) calloc (lc_MPOINTS, sizeof (double));
            for (iflux = 0; iflux < lc_MPOINTS; iflux++)
                lc_t_fit[iflux] = log10 (t_raw[iflux]);

            /* Verify monotonicity of the surviving knots */
            for (m = 1; m < lc_MPOINTS; m++) {
                if (lc_t_fit[m] <= lc_t_fit[m-1]) {
                    print ("  > WARNING: t_fit not monotone at valid knot %d: "
                           "log10(t)=%e <= %e\n", m, lc_t_fit[m], lc_t_fit[m-1]);
                }
            }

            /* Allocate M_UV_fit with the valid (possibly reduced) MPOINTS */
            lc_M_UV_fit = ARRAY_4D (lc_MPOINTS, NX3_TOT, NX2_TOT, NX1_TOT, double);
            lc_M_UV_all_zero = ARRAY_3D (NX3_TOT, NX2_TOT, NX1_TOT, int);

            icount = match = 0;

            while (fscanf (fptr, "%ld ", &ii) != EOF) {
                if (fscanf (fptr, "%ld %le %le", &jj, &x1in, &x2in) != 3) {
                    print ("! read_sirocco_fluxes: M_UV data read error\n");
                    QUIT_PLUTO(1);
                }

                /* Read all n_raw M values for this spatial cell into a temp buffer,
                 * then store only those whose t_fit column was valid.
                 * This keeps file parsing in sync even when rows are skipped.    */
                {
                    double *mbuf = (double *) malloc (n_raw * sizeof (double));
                    int found_cell = 0;

                    DOM_LOOP (k, j, i) {
                        if (fabs (1.0 - x1in / UNIT_LENGTH / x1[i]) < tol &&
                            fabs (1.0 - x2in / x2[j])              < tol) {
                            found_cell = 1;
                            for (iflux = 0; iflux < n_raw; iflux++) {
                                if (fscanf (fptr, "%le", &temp) != 1) {
                                    print ("! read_sirocco_fluxes: M_UV value read error\n");
                                    QUIT_PLUTO(1);
                                }
                                mbuf[iflux] = temp;
                            }
                            {
                                int n_floored = 0;  /* count bad knots for this cell */
                                int n_positive = 0;
                                for (iflux = 0; iflux < lc_MPOINTS; iflux++) {
                                    temp = mbuf[valid_idx[iflux]];
                                    /* Floor only genuinely unphysical M values (zero,
                                     * negative, non-finite) at 1e-4.  This avoids
                                     * log10(0) blow-up while keeping the floor well
                                     * below the electron-scattering limit (M=1), so
                                     * legitimate small-M cells from the CAK table are
                                     * left untouched.                                 */
                                    if (!isfinite (temp) || temp <= 0.0) {
                                        temp = 1.e-4;
                                        n_floored++;
                                    } else {
                                        n_positive++;
                                    }
                                    lc_M_UV_fit[iflux][k][j][i] = log10 (temp);
                                }
                                /* log10(M) cannot represent an exact zero.  Keep a
                                 * separate mask so an entirely zero CAK curve remains
                                 * physically M=0 in LineForce; the 1e-4 values above
                                 * are then only safe internal placeholders. */
                                lc_M_UV_all_zero[k][j][i] = (n_positive == 0);

                                if (lc_M_UV_all_zero[k][j][i]) {
                                    print ("  > NOTICE: M_UV cell (%ld,%ld): "
                                           "all %d knots non-positive; using M=0 "
                                           "(inspect CAK radiation/ion inputs)\n",
                                           ii, jj, lc_MPOINTS);
                                } else if (n_floored > 0) {
                                    print ("  > WARNING: M_UV cell (%ld,%ld): "
                                           "%d/%d knots non-positive; only those "
                                           "knots floored to 1e-4 for log interpolation\n",
                                           ii, jj, n_floored, lc_MPOINTS);
                                }
                            }
                            match = 1;
                            icount++;
                        }
                    }

                    if (!match) {
                        for (iflux = 0; iflux < n_raw; iflux++) {
                            if (fscanf (fptr, "%le", &mbuf[iflux]) != 1) {
                                print ("! read_sirocco_fluxes: M_UV skip error\n");
                                QUIT_PLUTO(1);
                            }
                        }
                    }
                    free (mbuf);
                }
                match = 0;
            }

            free (t_raw);
            free (valid_idx);
        }

        if (lc_MPOINTS > 0) {
            fclose (fptr);
            /* Every active local cell needs a mapped M(t) table. A unit or
               grid mismatch otherwise leaves uninitialized force values. */
            {
                int expected_cells = 0;
                DOM_LOOP(k, j, i) { expected_cells++; }
                if (icount != expected_cells) {
                    print("! read_sirocco_fluxes: M_UV_data.dat maps %d cells; expected %d on this rank. Check grid and coordinate units (r in cm, theta in radians).\n",
                          icount, expected_cells);
                    QUIT_PLUTO(1);
                }
            }
            {
                int icount_global = icount;
#ifdef PARALLEL
                MPI_Allreduce (&icount, &icount_global, 1, MPI_INT, MPI_SUM,
                               MPI_COMM_WORLD);
#endif
                print ("  > read %d M_UV cells (this proc; %d total)\n",
                       icount, icount_global);
            }
        }

        end_muv_block: ;   /* target for early-exit when no valid knots */

        } /* end if fptr != NULL */
    } /* end M_UV block */
}


/* ====================================================================
   2.  read_sirocco_heatcool
       flag=0 → first-call initialisation
       flag=1 → read sirocco_heatcool.dat and prefactors.dat
   ==================================================================== */
void read_sirocco_heatcool (Data *d, Grid *grid, int flag)
{
    int    i, j, k, ii, jj, nwords, icount;
    FILE   *fptr, *fptr_hc;
    char   aline[LC_LINELENGTH];
    double tol = 1e-5;
    double dens, rcen, thetacen;
    double comp_h_pre, comp_c_pre, xray_h_pre, brem_c_pre,
           line_c_pre, xi_ion_pre;
    double vol, t_e, t_r, xi, ne;
    double heat_xray, heat_comp, heat_lines, heat_ff;
    double cool_comp, cool_lines, cool_ff, n_h;
    double rho, r, nH, lx;

    /* ---- Initialisation pass (flag == 0) ---- */
    if (flag == 0) {

        /* Pre-populate sirocco_xi and sirocco_t_r with the analytical estimate
         * so that data.0000.dbl (written before any SIROCCO output exists) carries
         * a physically meaningful XI value.  pre_calc(0) in the Python driver
         * calls pload(0) and reads D.XI from data.0000.dbl; if xi were -1 here
         * the XI uservar would be written as -1 and pre_calc ratios would be
         * nonsensical.  All prefactors are set to unity (no correction yet).
         *
         * The analytical xi is also the correct fallback used by cooling.c and
         * userdef_output.c whenever lc_NFLUX_ANGLES == 0. */
        lx = g_inputParam[L_STAR] * g_inputParam[F_X];

        DOM_LOOP (k, j, i) {
#if COOLING != NO
            rho = d->Vc[RHO][k][j][i] * UNIT_DENSITY;
            r   = grid->x[IDIR][i]    * UNIT_LENGTH;
            nH  = rho / (1.43 * CONST_mp);

            d->sirocco_xi [k][j][i] = lx / (nH * r * r);  /* analytical estimate */
            d->sirocco_t_r[k][j][i] = g_inputParam[T_X];

            d->comp_h_pre[k][j][i] = 1.0;
            d->comp_c_pre[k][j][i] = 1.0;
            d->xray_h_pre[k][j][i] = 1.0;
            d->line_c_pre[k][j][i] = 1.0;
            d->brem_c_pre[k][j][i] = 1.0;
            d->xi_ion_pre[k][j][i] = 1.0;
#endif
        }
        print ("  > read_sirocco_heatcool: initialised (flag=0, analytical xi set)\n");
        return;
    }

    /* ---- Runtime pass (flag != 0): read sirocco_heatcool.dat ---- */
    fptr_hc = fopen ("sirocco_heatcool.dat", "r");

    if (fptr_hc == NULL) {
        print ("! read_sirocco_heatcool: required sirocco_heatcool.dat is missing\n");
        print ("! Run rad_hydro_files from the current SIROCCO memory2609 build before restarting PLUTO.\n");
        QUIT_PLUTO(1);
    } else {
        /* skip header */
        if (fgets (aline, LC_LINELENGTH, fptr_hc) == NULL) {
            print ("! read_sirocco_heatcool: sirocco_heatcool header read error\n");
            QUIT_PLUTO(1);
        }

        icount = 0;
        while (fgets (aline, LC_LINELENGTH, fptr_hc) != NULL) {
            nwords = sscanf (aline,
                "%d %d %le %le %le %le %le %le %le "
                "%le %le %le %le %le %le %le %le %le",
                &ii, &jj, &rcen, &thetacen, &vol, &t_e, &t_r, &xi, &ne,
                &heat_xray, &heat_comp, &heat_lines, &heat_ff,
                &cool_comp, &cool_lines, &cool_ff, &dens, &n_h);

            if (nwords != 18) {
                print ("! read_sirocco_heatcool: sirocco_heatcool format error "
                       "(nwords=%d)\n", nwords);
                QUIT_PLUTO(1);
            }

            DOM_LOOP (k, j, i) {
                if (fabs ((rcen / UNIT_LENGTH - grid->x[IDIR][i]) /
                          (rcen / UNIT_LENGTH)) < tol &&
                    fabs ((thetacen - grid->x[JDIR][j]) / thetacen) < tol) {

                    icount++;
#if COOLING != NO
                    d->sirocco_xi [k][j][i] = (xi  < 1.0)   ? 1.0   : xi;
                    d->sirocco_t_r[k][j][i] = (t_r < 1.e3)  ? 1.e3  : t_r;
#endif
                }
            }
        }
        fclose (fptr_hc);
        {
            int icount_global = icount;
#ifdef PARALLEL
            MPI_Allreduce (&icount, &icount_global, 1, MPI_INT, MPI_SUM,
                           MPI_COMM_WORLD);
#endif
            print ("  > read_sirocco_heatcool: read %d heatcool entries "
                   "(this proc; %d total)\n", icount, icount_global);
        }
    }

    /* ---- Read prefactors.dat ---- */
    fptr = fopen ("prefactors.dat", "r");

    if (fptr == NULL) {
        print ("  > read_sirocco_heatcool: no prefactors.dat, using unity\n");
        DOM_LOOP (k, j, i) {
#if COOLING != NO
            d->comp_h_pre[k][j][i] = 1.0;
            d->comp_c_pre[k][j][i] = 1.0;
            d->xray_h_pre[k][j][i] = 1.0;
            d->line_c_pre[k][j][i] = 1.0;
            d->brem_c_pre[k][j][i] = 1.0;
            d->xi_ion_pre[k][j][i] = 1.0;
#endif
        }
        return;
    }

    /* initialise to unity so missing cells keep sane values */
    DOM_LOOP (k, j, i) {
#if COOLING != NO
        d->comp_h_pre[k][j][i] = 1.0;
        d->comp_c_pre[k][j][i] = 1.0;
        d->xray_h_pre[k][j][i] = 1.0;
        d->line_c_pre[k][j][i] = 1.0;
        d->brem_c_pre[k][j][i] = 1.0;
        d->xi_ion_pre[k][j][i] = 1.0;
#endif
    }

    if (fgets (aline, LC_LINELENGTH, fptr) == NULL) {
        print ("! read_sirocco_heatcool: prefactors header read error\n");
        QUIT_PLUTO(1);
    }

    icount = 0;
    while (fgets (aline, LC_LINELENGTH, fptr) != NULL) {
        nwords = sscanf (aline,
            "%d %le %d %le %le %le %le %le %le %le %le",
            &ii, &rcen, &jj, &thetacen, &dens,
            &comp_h_pre, &comp_c_pre, &xray_h_pre,
            &brem_c_pre, &line_c_pre, &xi_ion_pre);

        if (nwords != 11) {
            print ("! read_sirocco_heatcool: prefactors format error "
                   "(nwords=%d)\n", nwords);
            QUIT_PLUTO(1);
        }

        DOM_LOOP (k, j, i) {
            if (fabs ((rcen / UNIT_LENGTH - grid->x[IDIR][i]) /
                      (rcen / UNIT_LENGTH)) < tol &&
                fabs ((thetacen - grid->x[JDIR][j]) / thetacen) < tol) {

                icount++;
#if COOLING != NO
                d->comp_h_pre[k][j][i] = comp_h_pre;
                d->comp_c_pre[k][j][i] = comp_c_pre;
                d->xray_h_pre[k][j][i] = xray_h_pre;
                d->line_c_pre[k][j][i] = line_c_pre;
                d->brem_c_pre[k][j][i] = brem_c_pre;
                d->xi_ion_pre[k][j][i] = xi_ion_pre;
#endif
            }
        }
    }

    fclose (fptr);
    {
        int icount_global = icount;
#ifdef PARALLEL
        MPI_Allreduce (&icount, &icount_global, 1, MPI_INT, MPI_SUM,
                       MPI_COMM_WORLD);
#endif
        print ("  > read_sirocco_heatcool: read %d prefactor entries "
               "(this proc; %d total)\n", icount, icount_global);
    }
}


/* ====================================================================
   3.  VGradCalc
       Computes |dv/ds| for each angular bin in each cell.
       Geometry arrays (lc_x11, lc_x22, offsets) are built only
       on the first call (lc_dvds_ready == 0).
   ==================================================================== */
#if LINE_DRIVEN_WIND != NO
void VGradCalc (const Data *d, Grid *grid)
{
    int    i, j, k, iangle;
    double v11[2], v12[2], v22[2], v21[2];
    double v1, v2, dx1, dx2;
    double x11[2], x22[2];
    double loc[2], ans1[2], ans2[2];
    double mod_flux, ds, maxds;
    double theta_angle;
    double x, z, vx1, vz1, vx2, vz2;

    double *x1  = grid->x [IDIR];
    double *x2  = grid->x [JDIR];

    /* ---- guard: no flux data loaded yet (cycle 0 not complete) ---- */
    if (lc_NFLUX_ANGLES == 0) {
        /* lc_dvds_array is NULL; LineForce will return zero — nothing to do */
        return;
    }

    /* ---- first-call geometry initialisation ---- */
    if (lc_dvds_ready == 0) {

        print ("  > VGradCalc: allocating geometry arrays\n");

        /* interpolation corner positions */
        lc_x11 = ARRAY_4D (2, NX3_TOT, NX2_TOT, NX1_TOT, double);
        lc_x22 = ARRAY_4D (2, NX3_TOT, NX2_TOT, NX1_TOT, double);

        /* dv/ds workspace */
        lc_dvds_r_off   = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);
        lc_dvds_t_off   = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);
        lc_dvds_mod_off = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);
        lc_dvds_array   = ARRAY_4D (lc_NFLUX_ANGLES, NX3_TOT, NX2_TOT, NX1_TOT, double);

        for (iangle = 0; iangle < lc_NFLUX_ANGLES; iangle++) {
            DOM_LOOP (k, j, i) {

                /* cell-face midpoints in CGS */
                lc_x11[0][k][j][i] = 0.5*(x1[i-1]+x1[i]) * UNIT_LENGTH;
                lc_x11[1][k][j][i] = 0.5*(x2[j-1]+x2[j]);
                lc_x22[0][k][j][i] = 0.5*(x1[i+1]+x1[i]) * UNIT_LENGTH;
                lc_x22[1][k][j][i] = 0.5*(x2[j+1]+x2[j]);

                /* path length to step for dv/ds estimate */
                maxds = fabs (lc_x22[0][k][j][i] - lc_x11[0][k][j][i]);
#if GEOMETRY == SPHERICAL
                {
                    double arc = x1[i]*UNIT_LENGTH *
                                 fabs (lc_x22[1][k][j][i] - lc_x11[1][k][j][i]);
                    if (arc < maxds) maxds = arc;
                }
#endif
                maxds *= 0.5;

                mod_flux = sqrt (  lc_flux_r[iangle][k][j][i]
                                 * lc_flux_r[iangle][k][j][i]
                                 + lc_flux_t[iangle][k][j][i]
                                 * lc_flux_t[iangle][k][j][i]
                                 + lc_flux_p[iangle][k][j][i]
                                 * lc_flux_p[iangle][k][j][i]);

                if (mod_flux == 0.0) {
                    lc_dvds_mod_off[iangle][k][j][i] = -999.0;
                } else {
#if GEOMETRY == SPHERICAL
                    x = x1[i]*sin(x2[j]) * UNIT_LENGTH;
                    z = x1[i]*cos(x2[j]) * UNIT_LENGTH;

                    theta_angle = (iangle + 0.5) * (2.0*CONST_PI) / 36.0;
                    dx1 = maxds * sin (theta_angle);
                    dx2 = maxds * cos (theta_angle);
                    ds  = sqrt (dx1*dx1 + dx2*dx2);

                    lc_dvds_r_off[iangle][k][j][i]   =
                        sqrt ((x+dx1)*(x+dx1) + (z+dx2)*(z+dx2));
                    lc_dvds_t_off[iangle][k][j][i]   =
                        atan ((x+dx1) / (z+dx2));
                    lc_dvds_mod_off[iangle][k][j][i] = ds;
#endif
                }
            }
        }

        print ("  > VGradCalc: geometry arrays ready\n");
        lc_dvds_ready = 1;
    }

    /* ---- compute dv/ds every time step ---- */
    for (iangle = 0; iangle < lc_NFLUX_ANGLES; iangle++) {
        DOM_LOOP (k, j, i) {

            /* vertex-averaged velocities */
            v11[0] = (d->Vc[VX1][k][j-1][i-1]+d->Vc[VX1][k][j-1][i]
                    + d->Vc[VX1][k][j  ][i-1]+d->Vc[VX1][k][j  ][i]) * 0.25;
            v11[1] = (d->Vc[VX2][k][j-1][i-1]+d->Vc[VX2][k][j-1][i]
                    + d->Vc[VX2][k][j  ][i-1]+d->Vc[VX2][k][j  ][i]) * 0.25;

            v12[0] = (d->Vc[VX1][k][j+1][i-1]+d->Vc[VX1][k][j  ][i-1]
                    + d->Vc[VX1][k][j+1][i  ]+d->Vc[VX1][k][j  ][i  ]) * 0.25;
            v12[1] = (d->Vc[VX2][k][j+1][i-1]+d->Vc[VX2][k][j  ][i-1]
                    + d->Vc[VX2][k][j+1][i  ]+d->Vc[VX2][k][j  ][i  ]) * 0.25;

            v22[0] = (d->Vc[VX1][k][j+1][i  ]+d->Vc[VX1][k][j+1][i+1]
                    + d->Vc[VX1][k][j  ][i+1]+d->Vc[VX1][k][j  ][i  ]) * 0.25;
            v22[1] = (d->Vc[VX2][k][j+1][i  ]+d->Vc[VX2][k][j+1][i+1]
                    + d->Vc[VX2][k][j  ][i+1]+d->Vc[VX2][k][j  ][i  ]) * 0.25;

            v21[0] = (d->Vc[VX1][k][j  ][i+1]+d->Vc[VX1][k][j-1][i+1]
                    + d->Vc[VX1][k][j-1][i  ]+d->Vc[VX1][k][j  ][i  ]) * 0.25;
            v21[1] = (d->Vc[VX2][k][j  ][i+1]+d->Vc[VX2][k][j-1][i+1]
                    + d->Vc[VX2][k][j-1][i  ]+d->Vc[VX2][k][j  ][i  ]) * 0.25;

            /* retrieve cached corners */
            x11[0] = lc_x11[0][k][j][i];
            x11[1] = lc_x11[1][k][j][i];
            x22[0] = lc_x22[0][k][j][i];
            x22[1] = lc_x22[1][k][j][i];

            /* velocity at cell centre */
            loc[0] = x1[i] * UNIT_LENGTH;
            loc[1] = x2[j];
            lc_bilinear (x11, x22, v11, v12, v21, v22, loc, ans1, 0);

            vx1 = (ans1[0]*UNIT_VELOCITY*sin(x2[j])
                +  ans1[1]*UNIT_VELOCITY*cos(x2[j]));
            vz1 = (ans1[0]*UNIT_VELOCITY*cos(x2[j])
                -  ans1[1]*UNIT_VELOCITY*sin(x2[j]));

            /* velocity at offset point */
            loc[0] = lc_dvds_r_off  [iangle][k][j][i];
            loc[1] = lc_dvds_t_off  [iangle][k][j][i];
            ds     = lc_dvds_mod_off[iangle][k][j][i];

            if (ds < 0.0) {                       /* -999 flag: no flux */
                lc_dvds_array[iangle][k][j][i] = -999.0;
            } else {
                lc_bilinear (x11, x22, v11, v12, v21, v22, loc, ans2, 0);

                vx2 = (ans2[0]*UNIT_VELOCITY*sin(loc[1])
                    +  ans2[1]*UNIT_VELOCITY*cos(loc[1]));
                vz2 = (ans2[0]*UNIT_VELOCITY*cos(loc[1])
                    -  ans2[1]*UNIT_VELOCITY*sin(loc[1]));

                theta_angle = (iangle + 0.5) * (2.0*CONST_PI) / 36.0;
                v1 = sin(theta_angle)*vx1 + cos(theta_angle)*vz1;
                v2 = sin(theta_angle)*vx2 + cos(theta_angle)*vz2;

                lc_dvds_array[iangle][k][j][i] = fabs ((v2-v1) / ds);
            }

            if (lc_dvds_array[iangle][k][j][i] != lc_dvds_array[iangle][k][j][i]) {
                print ("! VGradCalc: NaN at step=%ld iang=%d k=%d j=%d i=%d\n",
                       g_stepNumber, iangle, k, j, i);
                QUIT_PLUTO(1);
            }
        }
    }
}
#endif /* LINE_DRIVEN_WIND != NO */


/* ====================================================================
   4.  LineForce
       Computes radiation acceleration (r, theta, phi) at one cell.
   ==================================================================== */
void LineForce (double *v, double *grad,
                double x1, double x2, double x3,
                int i, int j, int k)
{
    int    iangle, ii;
    double *M_UV_arr = NULL;
    double sigma_e, rho, T, mu;
    double flux_r, flux_t, t_UV, M_UV, v_th;
    double M_max = 4400.0;

    /* Fixed sigma_T/m_p: matches cak_v4 and SIROCCO memory2609. */
    sigma_e  = 0.66524e-24 / 1.672661e-24;
    rho      = v[RHO] * UNIT_DENSITY;
    mu       = g_inputParam[MU];

#if EOS != ISOTHERMAL
    T = v[PRS] / v[RHO] * KELVIN * mu;
#else
    T = g_inputParam[T_ISO];
#endif

    /* copy per-cell M table if available (lc_MPOINTS > 0 after CAK run) */
    if (lc_MPOINTS > 0 && lc_M_UV_fit != NULL &&
        (lc_M_UV_all_zero == NULL || !lc_M_UV_all_zero[k][j][i])) {
        M_UV_arr = (double *) malloc (lc_MPOINTS * sizeof (double));
        for (ii = 0; ii < lc_MPOINTS; ii++)
            M_UV_arr[ii] = lc_M_UV_fit[ii][k][j][i];
    }

    v_th = sqrt (2.0 * CONST_kB * T / CONST_mp);

    grad[IDIR] = 0.0;
    grad[JDIR] = 0.0;
    grad[KDIR] = 0.0;

    for (iangle = 0; iangle < lc_NFLUX_ANGLES; iangle++) {

        flux_r = lc_flux_r[iangle][k][j][i];
        flux_t = lc_flux_t[iangle][k][j][i];

        if (lc_dvds_array[iangle][k][j][i] > 0.0) {
            t_UV = sigma_e * rho * v_th / lc_dvds_array[iangle][k][j][i];

            /* use tabulated M(t) if CAK has run, else pure e- scattering */
            if (M_UV_arr != NULL) {
                M_UV = lc_linterp (log10 (t_UV), lc_t_fit, M_UV_arr, lc_MPOINTS);
            } else {
                M_UV = 0.0;  /* no CAK table yet: only electron scattering */
            }
            if (M_UV > M_max) M_UV = M_max;
        } else {
            M_UV = 0.0;
        }

        grad[IDIR] += (1.0 + M_UV) * sigma_e * flux_r / CONST_c / UNIT_ACCELERATION;
        grad[JDIR] += (1.0 + M_UV) * sigma_e * flux_t / CONST_c / UNIT_ACCELERATION;
    }

    if (M_UV_arr) { free (M_UV_arr); M_UV_arr = NULL; }
}


/* ====================================================================
   Private helpers
   ==================================================================== */

/* ---- log-linear interpolation with boundary clamping ---- */
static double lc_linterp (double x, double *xarr, double *yarr, int n)
{
    int    idx = 0;
    double slope, expo;

    if (n <= 0) return 0.0;

    while (idx < n && xarr[idx] < x) idx++;

    if (idx == 0)   return pow (10.0, yarr[0]);
    if (idx >= n)   return pow (10.0, yarr[n-1]);

    if (!isfinite (xarr[idx-1]) || !isfinite (xarr[idx]) ||
        !isfinite (yarr[idx-1]) || !isfinite (yarr[idx])  ||
        xarr[idx] == xarr[idx-1]) {
        return 0.0;
    }

    slope = (yarr[idx] - yarr[idx-1]) / (xarr[idx] - xarr[idx-1]);
    expo  = yarr[idx-1] + slope * (x - xarr[idx-1]);
    return (isfinite (expo)) ? pow (10.0, expo) : 0.0;
}

/* ---- bilinear interpolation of 2-component velocity ---- */
static void lc_bilinear (double x11[2], double x22[2],
                          double v11[2], double v12[2],
                          double v21[2], double v22[2],
                          double test[2], double ans[2], int flag)
{
    double f1, f2, t1, t2;

    f1 = (test[0] - x11[0]) / (x22[0] - x11[0]);
    f2 = (test[1] - x11[1]) / (x22[1] - x11[1]);

    /* VX1 component */
    t1     = (1.0-f1)*v11[0] + f1*v21[0];
    t2     = (1.0-f1)*v12[0] + f1*v22[0];
    ans[0] = (1.0-f2)*t1 + f2*t2;

    /* VX2 component */
    t1     = (1.0-f1)*v11[1] + f1*v21[1];
    t2     = (1.0-f1)*v12[1] + f1*v22[1];
    ans[1] = (1.0-f2)*t1 + f2*t2;

    if (flag) {
        print ("  lc_bilinear: f1=%e f2=%e ans=(%e,%e)\n",
               f1, f2, ans[0], ans[1]);
    }
}

/* ===================================================================
   Accessor functions — expose lc_ static arrays to userdef_output.c
   without polluting the global namespace via pluto.h.
   Declare these in userdef_output.c with 'extern' before use.
   =================================================================== */

int    LC_GetNFluxAngles (void)                             { return lc_NFLUX_ANGLES; }
int    LC_GetMPoints     (void)                             { return lc_MPOINTS; }
int    LC_IsMUVAllZero   (int k, int j, int i)              { return lc_M_UV_all_zero ? lc_M_UV_all_zero[k][j][i] : 0; }
double LC_GetFluxR       (int a, int k, int j, int i)       { return lc_flux_r    ? lc_flux_r   [a][k][j][i] : 0.0; }
double LC_GetFluxT       (int a, int k, int j, int i)       { return lc_flux_t    ? lc_flux_t   [a][k][j][i] : 0.0; }
double LC_GetDvDs        (int a, int k, int j, int i)       { return lc_dvds_array? lc_dvds_array[a][k][j][i] : 0.0; }
double LC_GetMUVFit      (int m, int k, int j, int i)       { return lc_M_UV_fit  ? lc_M_UV_fit [m][k][j][i] : 0.0; }
double LC_GetTFit        (int m)                            { return lc_t_fit     ? lc_t_fit    [m]           : 0.0; }
double LC_LInterp        (double x, int n) {
    /* Wrapper around lc_linterp for use in userdef_output.c.
       Caller must build its own M_UV_arr from LC_GetMUVFit first. */
    (void)x; (void)n; return 0.0;  /* placeholder — caller uses own copy */
}
