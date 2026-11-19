/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Problem initialisation for the CV ideal disc-wind (cv_ideal).

  Sets up an axisymmetric 2-D (r, theta) Keplerian disc initial
  condition following Proga, Stone & Drew (1998).  Blondin
  cooling and SIROCCO line driving are enabled.

  The setup flag for VGradCalc geometry is now internal to
  line_connect.c (static lc_dvds_ready); no extern is needed here.

  \author  N. Higginbottom  (original)
  \author  A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date    April 2026
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/* ********************************************************************* */
{
    double cent_mass, rho_0, rho_alpha, disk_mdot, mu, dfloor;
    double rmin, gm_cgs, gm_code;
    double temp, cs, r_0, teff;
    double s, c, R, r, z;
    double v_k, R_1, tan_1;
    double rho_d, rho_a;


    /* ---- parameters from pluto.ini ---- */
    cent_mass = g_inputParam[CENT_MASS];
    disk_mdot = g_inputParam[DISK_MDOT];
    rho_0     = g_inputParam[RHO_0];
    rho_alpha = g_inputParam[RHO_ALPHA];
    mu        = g_inputParam[MU];
    dfloor    = g_inputParam[DFLOOR];

    rmin   = g_domBeg[IDIR];
    r_0    = rmin * UNIT_LENGTH;          /* inner radius in cgs */

    gm_cgs  = CONST_G * cent_mass;
    gm_code = gm_cgs / (UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY);

    r = x1;
    s = sin (x2);
    c = cos (x2);
    R = r * s;                            /* cylindrical radius */
    z = r * c;

    R_1   = (fabs (R)       > 1.e-12) ? 1.0 / R       : 0.0;
    tan_1 = (fabs (tan(x2)) > 1.e-12) ? 1.0 / tan(x2) : 0.0;

    /* effective disc temperature and sound speed */
    teff  = pow (3.0*gm_cgs*disk_mdot / (8.0*CONST_PI*CONST_sigma), 0.25);
    teff *= pow (r_0, -0.75);
    temp  = teff * pow (rmin / r, 0.75);

    cs  = sqrt (CONST_Rgas * temp / mu) / UNIT_VELOCITY;
    v_k = sqrt (gm_code / r);


    /* disc density profile (PSD98) */
    rho_0 = (rho_0 / UNIT_DENSITY) * pow (rmin / r, rho_alpha);
    rho_d = rho_0 * exp (-gm_code * tan_1 * tan_1 / (2.0 * cs * cs * r));
    rho_a = dfloor / UNIT_DENSITY;

    /* primitive variables */
    v[iVR]   = 0.0;
    v[iVTH]  = 0.0;
    v[iVPHI] = v_k * s;

    if (rho_d > rho_a) {
        v[RHO] = rho_d;
        v[TRC] = 1.0;
    } else {
        v[RHO] = rho_a;
        v[TRC] = 0.0;
    }

#if HAVE_ENERGY
    v[PRS]  = v[RHO] * temp / (KELVIN * mu);
    g_gamma = g_inputParam[GAMMA];
#endif

#if EOS == ISOTHERMAL
    temp = g_inputParam[T_ISO];
    g_isoSoundSpeed = sqrt (CONST_Rgas * temp / mu) / UNIT_VELOCITY;
#endif


    /* NOTE: lc_dvds_ready (the geometry-setup flag for VGradCalc) is a
       static variable inside line_connect.c.  It is initialised to 0
       by the C runtime and does NOT need to be reset here. */
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/* ********************************************************************* */
{
#if LINE_DRIVEN_WIND == SIROCCO_MODE && COOLING != NO
    /* Snapshot zero must contain finite unity thermal prefactors. */
    read_sirocco_heatcool(d,grid,0);
#endif
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* ********************************************************************* */
{
}


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 *  Assign user-defined boundary conditions.
 *
 *********************************************************************** */
{

	int   i, j, k, nv;
	double  *x1, *x2, *x3, *x2_glob;

	double r, theta, clight_code;
	double rho_0, rho_alpha, T_iso;
	double r_0, rho_mid, temp, mu;
	double dfloor;
	double cent_mass, gm_cgs, gm_code;
	double r_WD, teff, disk_mdot, rcyl;
	double theta_max;
	double tfloor, pfloor;
	double ss_fac;                 /* Shakura-Sunyaev inner-edge factor */


	x2_glob    = grid->x_glob[JDIR];

	rho_alpha  = g_inputParam[RHO_ALPHA];
	cent_mass  = g_inputParam[CENT_MASS];
  disk_mdot  = g_inputParam[DISK_MDOT];    /* disk accretion rate */
	T_iso      = g_inputParam[T_ISO];
	mu         = g_inputParam[MU];                      /* mean molecular weight */

	rho_0      = g_inputParam[RHO_0]  / UNIT_DENSITY;   /* density at r_0 - in the code unit */
	r_0        = g_inputParam[R_0]    / UNIT_LENGTH;    /* r_0 - typically set to R_IC in the code unit */
	dfloor     = g_inputParam[DFLOOR] / UNIT_DENSITY;   /* in code unit */

	tfloor = 5.e2;
	pfloor = dfloor * tfloor / (KELVIN * mu);


	r_WD      = g_domBeg[IDIR];
  theta_max = g_domEnd[JDIR];


	gm_cgs      = CONST_G * cent_mass;
	gm_code     = gm_cgs / (UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY);


	if (side == 0){

		/* gc is for conservative grid */
		x1 = grid->xgc[IDIR];
		x2 = grid->xgc[JDIR];
		x3 = grid->xgc[KDIR];

		/* avoid too small density near the boundary */
		RBox dom_box;
		TOT_LOOP(k,j,i){

			   int convert_to_cons = 0;


         /* ------------------------------------------------------------
            Identify the disc-midplane row by PROXIMITY to theta = pi/2,
            not by a hard-coded cell index.

            The previous test, "j == grid->np_int[JDIR] + 2", silently
            assumes nghost == 3 (SHOCK_FLATTENING MULTID forces this, so it
            happens to hold here) and that the domain stops exactly at
            theta = pi/2.  A proximity test (nearest cell centre to pi/2,
            within half a cell width) is correct regardless of nghost or
            where pi/2 falls in the domain.
            ------------------------------------------------------------ */
				 if (fabs(x2[j] - 0.5*CONST_PI) < 0.5*(x2[MIN(j+1,grid->np_tot[JDIR]-1)] - x2[MAX(j-1,0)]))  {

						r        = x1[i];  /*  code unit */
						theta    = x2[j];
						rcyl     = r * sin(theta);

						rho_mid  = rho_0 * pow((r/r_WD), -1.0 * rho_alpha);

						/* conserve momentum */
						d->Vc[VX2][k][j][i] = (d->Vc[RHO][k][j][i] * d->Vc[VX2][k][j][i]) / rho_mid;

						/* set the density at midplane - in the code unit */
						d->Vc[RHO][k][j][i] = rho_mid;

						/* set radial velocity at the midplane to zero */
						d->Vc[VX1][k][j][i] = 0.0;

						/* set rotatinal velocity at the midplane to the Keplarian */
						d->Vc[VX3][k][j][i] = sqrt(gm_code / r) * sin(theta);


						#if EOS != ISOTHERMAL
							teff  = pow(3.0 * gm_cgs * disk_mdot / (8.0 * CONST_PI * CONST_sigma), 0.25);
							teff *= pow(r_WD * UNIT_LENGTH, -0.75);

							/* (1 - sqrt(r_WD/rcyl)) is NEGATIVE in the radial ghost
							   zones (rcyl < r_WD) and pow(negative, 0.25) is NaN.
							   The NaN is overwritten by X1_BEG a moment later, so it
							   has never surfaced, but it is a live NaN in d->Vc[PRS]
							   in the meantime.  This guard touches ONLY ghost cells --
							   every active cell is bit-identical to the original. */
							ss_fac = 1.0 - sqrt(r_WD / rcyl);
							if (ss_fac < 0.0) ss_fac = 0.0;

							temp = teff * pow(r_WD / rcyl, 0.75) * pow(ss_fac, 0.25);

							d->Vc[PRS][k][j][i] = rho_mid * temp / (KELVIN * mu);
						#endif

						d->Vc[TRC][k][j][i] = 1.0;

				 }
                /* Apply the configured hydrodynamic density floor. */

				 if (d->Vc[RHO][k][j][i] < dfloor) {

					 d->Vc[RHO][k][j][i] = dfloor;

					 {
					   double v_k_loc = sqrt(gm_code / x1[i]);
					   double v_cap   = 2.0 * v_k_loc;

					   double vmag2 = d->Vc[VX1][k][j][i] * d->Vc[VX1][k][j][i]
					                + d->Vc[VX2][k][j][i] * d->Vc[VX2][k][j][i]
					                + d->Vc[VX3][k][j][i] * d->Vc[VX3][k][j][i];

					   if (vmag2 > v_cap * v_cap) {
					     double scale = v_cap / sqrt(vmag2);
					     d->Vc[VX1][k][j][i] *= scale;
					     d->Vc[VX2][k][j][i] *= scale;
					     d->Vc[VX3][k][j][i] *= scale;
					   }
					 }

					 #if EOS != ISOTHERMAL
					   temp = d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i] * KELVIN * mu;
					   if (temp < tfloor) {
					     d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * tfloor / (KELVIN * mu);
					   }
					 #endif

					 d->Vc[TRC][k][j][i] = 0.0;

					 convert_to_cons = 1;

				 }



			#if EOS != ISOTHERMAL

			 if (d->Vc[PRS][k][j][i] < pfloor) {
				 d->Vc[PRS][k][j][i] = pfloor;
				 convert_to_cons = 1;

			 }
 		  #endif

			 

			if (convert_to_cons) {
					 RBoxDefine (i, i, j, j, k, k, CENTER, &dom_box);
					 PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
				}




		  } /* DOM_LOOP() */
	  } /* if (side == 0) */


		if (side == X1_BEG){
	        /* gc is for conservative grid */
                x1 = grid->xgc[IDIR];
                x2 = grid->xgc[JDIR];


		  if (box->vpos == CENTER){
		    BOX_LOOP(box,k,j,i){
		      for (nv = 0; nv < NVAR; nv++){
		        d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
		      }
		      d->Vc[VX1][k][j][i] = MIN(d->Vc[VX1][k][j][i], 0.0);

		    }
		  }
		 }


		 if (side == X1_END){
			 if (box->vpos == CENTER){
				 BOX_LOOP(box,k,j,i){
					 for (nv = 0; nv < NVAR; nv++){
						 d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
					 }
					 d->Vc[VX1][k][j][i] = MAX(d->Vc[VX1][k][j][i], 0.0);
				 }
			 }
			}


			if (side == X2_BEG){                  /* -- X2_BEG (polar-axis) boundary -- */
				/* Keep the existing polar velocity reflection and scalar outflow. */
				if (box->vpos == CENTER){
					BOX_LOOP(box,k,j,i){
						for (nv = 0; nv < NVAR; nv++){
							d->Vc[nv][k][j][i] = d->Vc[nv][k][2*JBEG - j - 1][i];
						}

						d->Vc[VX2][k][j][i] *= -1.0;
						d->Vc[VX3][k][j][i] *= -1.0;   /* FIX: axisymmetric V_phi parity */


						d->Vc[RHO][k][j][i] = d->Vc[RHO][k][JBEG][i];         // apply outflow BCs
						#if EOS != ISOTHERMAL
								d->Vc[PRS][k][j][i] = d->Vc[PRS][k][JBEG][i];     // apply outflow BCs
						#endif
					}
				}
			}



}



#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector (double *v, double *g,
                      double x1, double x2, double x3)
/* ********************************************************************* */
{
    double cent_mass = g_inputParam[CENT_MASS];
    double gm_cgs    = CONST_G * cent_mass;
    double gm_code   = gm_cgs / (UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY);

    g[IDIR] = -gm_code / (x1 * x1);
    g[JDIR] =  0.0;
    g[KDIR] =  0.0;
}

/* ********************************************************************* */
double BodyForcePotential (double x1, double x2, double x3)
/* ********************************************************************* */
{
    return 0.0;
}
#endif
