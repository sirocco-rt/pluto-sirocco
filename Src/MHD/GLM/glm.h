/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Header file for GLM Divergence Cleaning

  Contains function prototypes and global variable declaration
  for the GLM formulation to control the divergence-free condition 
  of magnetic field.

  \authors A. Mignone (andrea.mignone@unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Sep 03, 2020

  \b References
     - "A Second-order unsplit Godunov scheme for cell-centered MHD:
        The CTU-GLM scheme"\n
        Mignone \& Tzeferacos, JCP (2010) 229, 2117

     - "High-order conservative finite difference GLM-MHD scheme for
        cell-centered MHD"\n
        Mignone, Tzeferacos \& Bodo, JCP (2010) 229, 5896
        
*/
/* ///////////////////////////////////////////////////////////////////// */
#define GLM_MHD

#ifndef GLM_ALPHA
  #define GLM_ALPHA    0.1  /**< Sets the damping rate of monopoles. */
#endif

#ifndef GLM_EXTENDED
 #define GLM_EXTENDED  NO  /**< The GLM_EXTENDED macro may be turned to YES to enable
                                the extended GLM formalism. 
                                Although it breaks conservation of momentum and
                                energy, it has proven to be more robust in treating
                                low-beta plasma. */
#endif

#ifndef GLM_COMPUTE_DIVB
  #define GLM_COMPUTE_DIVB  NO
#endif

#if (PHYSICS == ResRMHD) && !(defined GLM_COMPUTE_DIVE)
  #if CHARGE_SCHEME == 1
    #define GLM_COMPUTE_DIVE  YES
  #else
    #define GLM_COMPUTE_DIVE  NO
  #endif
#else
  #define GLM_COMPUTE_DIVE  NO
#endif

/* with chombo, COMPUTE_DIVB must be 
   disabled or a segfault will occur */

#ifdef CHOMBO
 #undef GLM_COMPUTE_DIVB  
 #undef GLM_COMPUTE_DIVE
 #define GLM_COMPUTE_DIVB NO
 #define GLM_COMPUTE_DIVE NO
#endif

extern double glm_ch; /**< The propagation speed of divergence error. */
    
void  GLM_Solve (const Sweep *, int, int, Grid *);
void  GLM_Init      (const Data *, const timeStep *, Grid *);
void  GLM_Source (const Data *, double, Grid *);
void  GLM_ExtendedSource (const Sweep *, double, int, int, Grid *);

#if GLM_COMPUTE_DIVB == YES
 void GLM_ComputeDivB(const Sweep *sweep, Grid *grid);
 double ***GLM_GetDivB(void);
#endif

#if GLM_COMPUTE_DIVE == YES
 void GLM_ComputeDivE(const Sweep *sweep, Grid *grid);
 double ***GLM_GetDivE(void);
#endif

