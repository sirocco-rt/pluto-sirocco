/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief This file (`init.c`) contains user-supplied functions for
         problem configuration and initialization of the Line-Driven
         Disk Wind simulation.

  \details The provided functions are used to set up and configure
           the initial conditions, boundary conditions, and other
           problem-specific parameters.

  \author Nick H
  \date   December 3, 2018

  \revised_by A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date       December 31, 2024
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"

extern double gFlux_r, gFlux_t, gFlux_p;

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{

	 double cent_mass, rho_alpha, disk_mdot;
 	 double rho_0, r_0, tx, T_iso, mu, dfloor;
 	 double r, temp, cs2_cgs, rho_cgs, r_cyl;
 	 double gm_cgs, rho_d, rho_a, prs_a, prs_d;
	 double temp_star, teff, r_WD, cs2_a;

   /* constant parameters from pluto.ini */

 	 cent_mass = g_inputParam[CENT_MASS];    /* central mass */
 	 rho_alpha = g_inputParam[RHO_ALPHA];    /* drop off exponent of density */
 	 disk_mdot = g_inputParam[DISK_MDOT];    /* disk accretion rate */
 	 rho_0     = g_inputParam[RHO_0];        /* density at r_0 set to R_ISCO */
   tx        = g_inputParam[T_x];          /* temperature of xray source  */
   T_iso     = g_inputParam[T_ISO];        /* isothermal temperature  */
 	 mu        = g_inputParam[MU];           /* mean molecular weight */
   dfloor    = g_inputParam[DFLOOR];       /* minimum density */

   r_0       = g_domBeg[IDIR] * UNIT_LENGTH;          /* r_0 - typically set to R_ISCO */

	 /* constants conversion for cgs units */
	 gm_cgs  = CONST_G * cent_mass;
	 r       = x1 * UNIT_LENGTH;

	 r_cyl   = r * sin(x2);
	 r_WD    = g_domBeg[IDIR];


   /* set gamma value for the equation of state */

   #if EOS == ISOTHERMAL
			temp = T_iso;
			g_isoSoundSpeed = sqrt(CONST_Rgas * temp / 0.6) / UNIT_VELOCITY;
			cs2_cgs = g_isoSoundSpeed * g_isoSoundSpeed * UNIT_VELOCITY * UNIT_VELOCITY;
   #elif EOS == IDEAL
	 		teff  = pow(3.0 * gm_cgs * disk_mdot / (8.0 * CONST_PI * CONST_sigma), 0.25);
			teff *= pow(r_WD * UNIT_LENGTH, -0.75);
//
			temp   = teff * pow(r_WD / x1, 0.75);
			cs2_cgs = (CONST_Rgas * temp / 0.6);
   #endif

  /* calculate density in cgs units - the expression from PSD98 */
  rho_cgs = rho_0 * exp(-1.0 * gm_cgs / (2.0 * cs2_cgs * r * tan(x2) * tan(x2)));


  rho_a = dfloor  / UNIT_DENSITY;
  rho_d = rho_cgs / UNIT_DENSITY;


  v[VX1] = v[VX2] = v[VX3] = 0.0;

	if (rho_d > rho_a)  {
 	  v[RHO] = rho_d;
    v[TRC] = 1.0;
  } else {
    v[RHO] = rho_a;
    v[TRC] = 0.0;
   }


 	v[VX3] = (sqrt((CONST_G * cent_mass ) / r) * sin(x2)) / UNIT_VELOCITY;


 	/* this converts temperature (in K) to pressure in code units */
 	#if HAVE_ENERGY
 		 v[PRS] = v[RHO] * temp / (KELVIN * mu);
 	#endif

  #if LINE_DRIVEN_WIND != NO
    dvds_setup_flag = 0;
  #endif



}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif


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
	double dfloor, dden, cs, dfact;
	double cent_mass, gm_cgs, gm_code;
	double r_WD, teff, disk_mdot, rcyl;
	double theta_max;
	double tfloor, pfloor;


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

				 // Check if density is negative or below the floor

				 if (d->Vc[RHO][k][j][i] < dfloor) {

					 double rho_new = (d->Vc[RHO][k][j][i] > 0) ? d->Vc[RHO][k][j][i] : 0.0;
					 if (rho_new < dfloor) {
							 rho_new = dfloor;  // Set density to floor if it’s negative or too low
					 }

					 if ( d->Vc[RHO][k][j][i] < 0.0 ) d->Vc[RHO][k][j][i] = dfloor;

					 #if EOS != ISOTHERMAL
						 cs    = sqrt(g_gamma * d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i]);
					 #endif


					 dfact = d->Vc[RHO][k][j][i] / dfloor;


					 d->Vc[RHO][k][j][i] = dfloor;

					 /* -----------------------------------------------------------------------------
					 ........... To conserve momentum and Energy, we modify velocities .............
					 // -------------------------------------------------------------------------- */

					 d->Vc[VX1][k][j][i] = dfact * d->Vc[VX1][k][j][i];
					 d->Vc[VX2][k][j][i] = dfact * d->Vc[VX2][k][j][i];
					 d->Vc[VX3][k][j][i] = dfact * d->Vc[VX3][k][j][i];


 //				 Ensure temperaure does not fall below the temperaure floor
					 #if EOS != ISOTHERMAL
//
					 	 d->Vc[PRS][k][j][i] = pow(cs, 2) * d->Vc[RHO][k][j][i] / g_gamma;

						 temp = d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i] * KELVIN * mu;

						 if (temp < tfloor) {
							 temp = tfloor;
						 	 d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * temp / (KELVIN * mu);
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


				 /* this should be the last 'real' theta bin - before the ghost zones. */

				 if (j == grid->np_int[JDIR] + 2 && 	(fabs(x2_glob[grid->np_int_glob[JDIR]+2]
								- x2[j]) / x2_glob[grid->np_int_glob[JDIR] + 2]) < 1e-6)  {

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
							temp = teff * pow(r_WD / rcyl, 0.75) * pow(1.0 - sqrt(r_WD / rcyl) ,0.25);
//
							d->Vc[PRS][k][j][i] = rho_mid * temp / (KELVIN * mu);
						#endif

						d->Vc[TRC][k][j][i] = 1.0;

				 }


		  } /* DOM_LOOP() */
	  } /* if (side == 0) */


		if (side == X1_BEG){
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


			if (side == X2_BEG){                  /* -- X2_BEG boundary --  Hybrid BC (Reflective for Velocity, Outflow for Scalars) */
				if (box->vpos == CENTER){
					BOX_LOOP(box,k,j,i){
						for (nv = 0; nv < NVAR; nv++){
							d->Vc[nv][k][j][i] = d->Vc[nv][k][2*JBEG - j - 1][i];
						}

						d->Vc[VX2][k][j][i] *= -1.0;

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
		void BodyForceVector(double *v, double *g, double x1, double x2, double x3)

		{

			double cent_mass, gm_cgs, gm_code;

			cent_mass = g_inputParam[CENT_MASS];
			gm_cgs    = CONST_G * cent_mass;
			gm_code   = gm_cgs / (UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY);


			g[IDIR] = - 1.0 * gm_code / (x1 * x1);
			g[JDIR] = 0.0;
			g[KDIR] = 0.0;

		}

		/* ********************************************************************* */
		double BodyForcePotential(double x1, double x2, double x3)
		/*!
		 * Return the gravitational potential as function of the coordinates.
		 *
		 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
		 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
		 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
		 *
		 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
		 *
		 *********************************************************************** */
		{
		}


#endif
