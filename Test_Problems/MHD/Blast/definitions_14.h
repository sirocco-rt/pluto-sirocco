#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  P_IN                           0
#define  P_OUT                          1
#define  BMAG                           2
#define  THETA                          3
#define  PHI                            4
#define  RADIUS                         5
#define  GAMMA                          6

/* [Beg] user-defined constants (do not change this line) */

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CHAR_LIMITING                  NO
#define  LIMITER                        VANLEER_LIM
#define  CT_EMF_AVERAGE                 UCT_HLLD
#define  CT_EN_CORRECTION               YES

/* [End] user-defined constants (do not change this line) */
