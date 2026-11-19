/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  PLUTO-SIROCCO line-driven wind connection module.

  This file is the canonical location of the subroutines that couple
  PLUTO to the SIROCCO Monte Carlo radiative-transfer code and the
  CAK force-multiplier tables.  It handles:

    - read_sirocco_fluxes()    : reads directional_flux_r/theta/phi.dat
                                 and M_UV_data.dat written by SIROCCO/CAK
    - read_sirocco_heatcool()  : reads py_heatcool.dat and prefactors.dat
    - VGradCalc()              : computes |dv/ds| for each angular bin
    - bilinear()               : bilinear interpolation helper
    - linterp()                : log-linear interpolation helper
    - LineForce()              : computes radiation acceleration vector

  The full implementation of these routines will be introduced in a
  later step.  For now this file provides only the skeleton so that
  the LineDriven/ directory and its makefile are correctly wired into
  the PLUTO build system.

  \author  A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date    April 2026
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ------------------------------------------------------------------ */
/*  Stub implementations — will be replaced in a later step           */
/* ------------------------------------------------------------------ */

void read_sirocco_fluxes (Data *d, Grid *grid)
{
  print ("  > read_sirocco_fluxes: not yet implemented\n");
}

void read_sirocco_heatcool (Data *d, Grid *grid, int flag)
{
  print ("  > read_sirocco_heatcool: not yet implemented\n");
}

#if LINE_DRIVEN_WIND != NO
void VGradCalc (const Data *d, Grid *grid)
{
  print ("  > VGradCalc: not yet implemented\n");
}
#endif

void LineForce (double *v, double *grad, double x1, double x2, double x3,
                int i, int j, int k)
{
  grad[IDIR] = 0.0;
  grad[JDIR] = 0.0;
  grad[KDIR] = 0.0;
}
