/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Header for the BLONDIN cooling/heating module.

  Declares the main entry point BlondinCooling() used by split_source.c,
  and the ne_rat() helper that returns the electron-to-hydrogen number
  density ratio as a function of ionization parameter xi.

  The full physics is documented in the PLUTO-SIROCCO implementation
  notes (Mosallanezhad et al. 2025, MNRAS 541, 2393).

  \author  A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date    April 2026
*/
/* ///////////////////////////////////////////////////////////////////// */
#ifndef BLONDIN_COOLING_H
#define BLONDIN_COOLING_H

void   BlondinCooling (Data_Arr, const Data *, double, timeStep *, Grid *);
double ne_rat          (double xi);

#endif /* BLONDIN_COOLING_H */
