#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            2

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  MACH                           0
#define  NU_VISC                        1

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MC_LIM
#define  CHOMBO_REF_VAR                 RHO

/* [End] user-defined constants (do not change this line) */
