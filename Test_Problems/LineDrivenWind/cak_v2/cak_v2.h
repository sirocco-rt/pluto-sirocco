#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>



#define nelements 30
#define log_t_min -10.0
#define log_t_max 7
#define n_t 200
#define LINELENGTH 2000


#define IONMODE_ML93 3          // Lucy Mazzali
#define IONMODE_MATRIX_BB 8     // matrix solver BB model
#define IONMODE_MATRIX_SPECTRALMODEL 9  // matrix solver spectral model based on power laws


#define ARRAY_2D(nx,ny,type)       (type   **)Array2D(nx,ny,sizeof(type))

//ion data - all of size nion

int *ion_info_z;          //atomic number
int *ion_info_state;      //state
int *ion_info_nlevels;    //number of levels
int *ion_info_nlines;     //number of lines
double *ion_info_xi;      //ionization state - not used  


//int ion_info_z[NIONS];          //atomic number
//int ion_info_state[NIONS];      //state
//int ion_info_nlevels[NIONS];    //number of levels
//int ion_info_nlines[NIONS];
//double ion_info_xi[NIONS];      //ionization state - not used  



//cell data

//double ion_fracs[nelements][nelements+1];
double *T_e, *nnh, *nne, *rho, *v_thermal, *sigma_e;
//double T_e[NCELLS], nnh[NCELLS], nne[NCELLS], rho[NCELLS], v_thermal[NCELLS], sigma_e[NCELLS];
int *cell_j; 
int *cell_i;



double **t;

double *param1; //hc/4pi/rho/vth


double T_max; //the maximum temperature seen in the model - used to compute exponential table
double T_min; //the minimum temperature seen in the model - used to compute exponential table
double line_fmin; //the minimum line frequency
double line_fmax; //the maximum line frequency

int nions,nline_tot;
double **ion_fracs;

double **by_ion;

//2D arrays for level data of size n_ions x n_levels

double **lev_energy;
double **lev_weight;
double **lev_pops;



int *trans_ion;
int *trans_lower;
int *trans_upper;
double *trans_freq;
double *trans_lfreq;
double *trans_A_ul;
double *trans_B_ul;
double *trans_B_lu;
int *trans_iband; //The band of the spectral model that this line lies in


/* Things to do with spectral model */

int model_type;
int nbands;

double **f1;
double **f2;
int **model;
double **pl_w;
double **pl_alpha;
double **exp_w;
double **exp_temp;

double *t_r;
double *W;

double *mod_fmin;
double *mod_fmax;
double *band_limits;

double *J;


#define SPEC_MOD_PL  1
#define SPEC_MOD_EXP  2
#define SPEC_MOD_FAIL -1





double **M_UV_array;

int ncells;

double g_usedMemory;

int npts_cont;
double *cont, *cont_nu;
double cont_norm;



double *level_energy;
double *level_pops;
double *level_weights;


int read_ionfs();
int read_cont();
int read_orig_cont();
double model_jnu(double freq, double lfreq,int icell,int iband,int debug);
double int_jnu(int icell, int iband);
double exp1(double x); 
double exp_lookup_init(); 
double num_int(double (*func)(double, void *), double a, double b, double eps);
double bb(double freq, void *params);
double bb2(double freq, double w, double t_r);

int read_line_data();
double get_ff(double nu_line);
char **Array2D (int, int, size_t);



#define MH 1.6733e-24 	
#define	KB 1.380658e-16
#define Ryd 2.1798741e-11
#define H 6.6260755e-27
#define C 2.99792458e10
#define EV 1.6021772e-12
#define CLIGHTSQUAREDOVERTWOH 6.7819570e46
#define HCLIGHTOVERFOURPI 1.580764662876770e-17
#define SIGMA_T 6.6524e-25 /* Thomson cross-section */
