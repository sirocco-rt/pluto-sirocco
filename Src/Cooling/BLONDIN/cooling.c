/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  BLONDIN cooling/heating module for PLUTO-SIROCCO.

  Implements a modified Blondin (1994) formula for radiative
  cooling/heating in X-ray ionised disc winds.  Five physical processes
  are included:

    1. Compton heating
    2. Compton cooling
    3. X-ray photoionisation heating
    4. Bremsstrahlung cooling
    5. Line cooling

  The ionisation parameter xi and radiation temperature T_x used here
  are supplied by SIROCCO via the read_sirocco_heatcool() routine in
  line_connect.c.  K-factor prefactors (comp_h_pre, comp_c_pre, etc.)
  are computed by pluto_sirocco_sub.py and stored in prefactors.dat.

  The full physics will be introduced in a later step.  This file
  currently provides the skeleton so that the BLONDIN/ directory and
  its makefile compile cleanly.

  \author  A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date    April 2026

  \b References
    - Blondin (1994), ApJ 435, 756
    - Mosallanezhad et al. (2025), MNRAS 541, 2393
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "cooling.h"

/* ------------------------------------------------------------------ */
/*!
 * \brief  Electron-to-hydrogen number density ratio.
 *
 * Returns ne/nH as a piecewise linear fit to photoionisation
 * equilibrium, as a function of log10(xi).  Will be fully
 * implemented in a later step.
 *
 * \param [in]  xi   Ionisation parameter (erg cm s^-1)
 * \return            ne/nH
 */
/* ------------------------------------------------------------------ */
double ne_rat (double xi)
{
  return 1.21;   /* placeholder: fully ionised solar abundance */
}

/* ------------------------------------------------------------------ */
/*!
 * \brief  Modified Blondin cooling/heating source term.
 *
 * Called by split_source.c in operator-split fashion.
 * Full implementation will be provided in a later step.
 */
/* ------------------------------------------------------------------ */
void BlondinCooling (Data_Arr VV, const Data *d, double dt,
                     timeStep *Dts, Grid *grid)
{
  print ("  > BlondinCooling: not yet implemented\n");
}
