#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            PVTE_LAW
#define  ENTROPY_SWITCH                 NO
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ETA                            0
#define  MACH                           1
#define  TJET                           2

/* [Beg] user-defined constants (do not change this line) */

#define  INCLUDE_JDIR                   NO
#define  LIMITER                        OSPRE_LIM
#define  PV_TEMPERATURE_TABLE           NO
#define  TV_ENERGY_TABLE                NO
#define  UNIT_DENSITY                   (1.e3*CONST_amu)
#define  UNIT_LENGTH                    2.5e15
#define  UNIT_VELOCITY                  1.e5
#define  WARNING_MESSAGES               NO

/* [End] user-defined constants (do not change this line) */
