#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            2

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      RK_LEGENDRE
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  OMEGA                          0
#define  REYN                           1

/* [Beg] user-defined constants (do not change this line) */

#define  INCLUDE_JDIR                   NO
#define  LIMITER                        MC_LIM

/* [End] user-defined constants (do not change this line) */
