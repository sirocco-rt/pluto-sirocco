#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENOZ
#define  TIME_STEPPING                  RK4
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define  RHO_IN                         0
#define  PRS_IN                         1
#define  RHO_OUT                        2
#define  PRS_OUT                        3
#define  BMAG                           4
#define  THETA                          5
#define  PHI                            6
#define  RADIUS                         7

/* [Beg] user-defined constants (do not change this line) */

#define  ASSIGN_VECTOR_POTENTIAL        NO
#define  CHAR_LIMITING                  NO
#define  CT_EMF_AVERAGE                 UCT_GFORCE
#define  CT_EN_CORRECTION               YES
#define  HO_LAP_LIMITER                 HO_RDER_LIM
#define  HO_LAP_LIMITER_BUF             0
#define  HO_LAP_LIMITER_EPS             0.2
#define  HO_ORDER_REDUCTION             LINEAR
#define  HO_PRIMITIVE_BC                YES
#define  HO_STARTUP_NGAUSS              4
#define  LIMITER                        VANLEER_LIM
#define  QUIT_ON_FIX                    NO
#define  SHOCK_FLATTENING               NO
#define  WARNING_MESSAGES               YES

/* [End] user-defined constants (do not change this line) */
