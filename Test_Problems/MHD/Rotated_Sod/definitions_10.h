#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENOZ
#define  TIME_STEPPING                  RK4
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            20

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

#define  RHO_LEFT                       0
#define  VX_LEFT                        1
#define  VY_LEFT                        2
#define  VZ_LEFT                        3
#define  BY_LEFT                        4
#define  BZ_LEFT                        5
#define  PR_LEFT                        6
#define  RHO_RIGHT                      7
#define  VX_RIGHT                       8
#define  VY_RIGHT                       9
#define  VZ_RIGHT                       10
#define  BY_RIGHT                       11
#define  BZ_RIGHT                       12
#define  PR_RIGHT                       13
#define  BX_CONST                       14
#define  GAMMA_EOS                      15
#define  INT_JX                         16
#define  INT_JY                         17
#define  INT_KX                         18
#define  INT_KZ                         19

/* [Beg] user-defined constants (do not change this line) */

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CHAR_LIMITING                  YES
#define  CHECK_DIVB_CONDITION           NO
#define  CT_EMF_AVERAGE                 UCT_HLLD
#define  CT_EN_CORRECTION               NO
#define  DIVIDE_BY_4PI                  YES
#define  HO_LAP_LIMITER                 HO_RDER_LIM
#define  HO_LAP_LIMITER_EPS             0.1
#define  HO_ORDER_REDUCTION             WENO3
#define  HO_PRIMITIVE_BC                NO
#define  HO_STARTUP_NGAUSS              4
#define  LIMITER                        VANLEER_LIM
#define  NGHOST_USR                     7
#define  SHOCK_FLATTENING               NO
#define  WARNING_MESSAGES               YES

/* [End] user-defined constants (do not change this line) */
