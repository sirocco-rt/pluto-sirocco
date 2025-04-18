/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Underexpanded jet.

  The domain is initialized with a static medium, whereas in the bottom
  boundary the following conditions are applied (see Section 5.2 of [Mig07])
  \f[
    P_{\rm jet}   = \frac{P_{\rm ratio}}{\Gamma}\left(\frac{2}{\Gamma + 1}\right)^{\Gamma/(\Gamma - 1)}
    \,,\quad
    \rho_{\rm jet} = \rho_{\rm ratio}\left(\frac{2}{\Gamma + 1}\right)^{1/(\Gamma - 1)}
    \,,\quad
    v_{\rm jet}   = \sqrt{\frac{\Gamma P_{\rm jet}}{\rho_{\rm jet}}}\,.
  \f]

  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[DN_RATIO]</tt>: sets the value of \f$\rho_{\rm ratio}\f$;
  - <tt>g_inputParam[PR_RATIO]</tt>: sets the value of \f$P_{\rm ratio}\f$;

  Configurations:

  - #01: Second order accuracy;
  - #02: Third order accuracy;

  \author A. Mignone (andrea.mignone@unito.it)
  \date   March 2, 2017

  \b References: 
     - [Mig07]: "PLUTO: a numerical code for computational astrophysics",
       Mignone et al., ApJS (2007), 170, 228

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  g_gamma = 5./3.;

  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = 1.0/g_gamma;
  us[TRC] = 0.0;

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
/* 
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  *R;
  double  pjet, dnjet, vjet;
  double  scrh;

  scrh = 1.0/(g_gamma - 1.0);

  R     = grid->xgc[IDIR];
  pjet  = g_inputParam[PR_RATIO]*pow(2.0/(g_gamma + 1.0),g_gamma*scrh)/g_gamma;
  dnjet = g_inputParam[DN_RATIO]*pow(2.0/(g_gamma + 1.0),scrh);
  vjet  = sqrt(g_gamma*pjet/dnjet);

#if GEOMETRY == CYLINDRICAL
  if (side == X2_BEG){
    BOX_LOOP(box,k,j,i){
      if (R[i] <= 1.) {
        d->Vc[RHO][k][j][i] = dnjet;
        d->Vc[VX1][k][j][i] = 0.;
        d->Vc[VX2][k][j][i] = vjet;
        d->Vc[PRS][k][j][i] = pjet;
      } else {
        d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][2*JBEG - j - 1][i];
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][2*JBEG - j - 1][i];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][2*JBEG - j - 1][i];
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][2*JBEG - j - 1][i];
      }
    }
  } 
#elif GEOMETRY == POLAR
  if (side == X3_BEG){
    BOX_LOOP(box,k,j,i){
      if (R[i] <= 1.) {
        d->Vc[RHO][k][j][i] = dnjet;
        d->Vc[VX1][k][j][i] = 0.;
        d->Vc[VX2][k][j][i] = 0.;
        d->Vc[VX3][k][j][i] = vjet;
        d->Vc[PRS][k][j][i] = pjet;
      } else {
        d->Vc[RHO][k][j][i] =  d->Vc[RHO][2*KBEG - k - 1][j][i];
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][2*KBEG - k - 1][j][i];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][2*KBEG - k - 1][j][i];
        d->Vc[VX3][k][j][i] = -d->Vc[VX3][2*KBEG - k - 1][j][i];
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][2*KBEG - k - 1][j][i];
      }
    }
  }
#endif
}

