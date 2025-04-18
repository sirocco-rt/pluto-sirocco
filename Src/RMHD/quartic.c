/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Solve quartic and cubic equations.



#####################################################################
  THIS FILE HAS BEEN DISMISSED ON October 4 2016.
  RELEVANT FUNCTIONS INCORPORATED INSIDE math_root_finders.c
#####################################################################

  

  \author A. Mignone (andrea.mignone@unito.it)
  \date   Oct 4, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void CheckSolutions(double b, double c, double d, double e, double *z);
void PrintSolution (double *z);
void QuarticPrintCoeffs(double b, double c, double d, double e);
int CubicSolve2 (double b, double c, double d, double z[]);

static int debug_print=0;

/* ********************************************************************* */
int QuarticSolve (double b, double c, double d, double e, double *z)
/*!
 * Solve a quartic equation in the form 
 * \f[
 *      z^4 + bz^3 + cz^2 + dz + e = 0
 * \f]
 *  For its purpose, it is assumed that \b ALL roots are double. 
 *  This makes things faster.
 *
 * \param [in] b   coefficient of the quartic
 * \param [in] c   coefficient of the quartic
 * \param [in] d   coefficient of the quartic
 * \param [in] e   coefficient of the quartic
 * \param [out] z  vector containing the 
 *                 (double) roots of the quartic
 *   
 * \b Reference:
 *
 *   http://www.1728.com/quartic2.htm 
 * 
 * \return Return 0 on success, 1 if cubic solver fails, 2 if NaN has
 *         been found and 3 if quartic is not satisfied.
 *
 *********************************************************************** */
{
  int status;
  status = QuarticSolveNew(b,c,d,e,z);

  if (status != 0){
    debug_print = 1;
    QuarticSolveNew(b,c,d,e,z);
    QUIT_PLUTO(1);
  }
}
  

/* ********************************************************************* */
int QuarticSolveNew (double b, double c, double d, double e, double *z)
/*!
 *
 *
 *********************************************************************** */
{
  int  n, j, ifail;
  double b2, sq;
  double a2, a1, a0, u[4];
  double p, q, r, f;
  const double three_256 = 3.0/256.0;
  const double one_64 = 1.0/64.0;
  double sz1, sz2, sz3, sz4;

  b2 = b*b;

/* --------------------------------------------------------------
   1) Compute cubic coefficients using the method outlined in
      http://eqworld.ipmnet.ru/En/solutions/ae/ae0108.pdf    
   -------------------------------------------------------------- */

  p = c - b2*0.375;
  q = d + b2*b*0.125 - b*c*0.5;
  r = e - 3.0*b2*b2/256.0 + b2*c/16.0 - 0.25*b*d;
  
  a2 = 2.0*p;
  a1 = p*p - 4.0*r;
  a0 = -q*q;

  ifail = CubicSolve(a2, a1, a0, u);
  if (ifail != 0) return 1;

  u[0] = MAX(u[0],0.0);
  u[1] = MAX(u[1],0.0);
  u[2] = MAX(u[2],0.0);


  if (u[0] != u[0] || u[1] != u[1] || u[2] != u[2]) return 1;

  sq  = -0.5*DSIGN(q); 
  sz1 = sq*sqrt(u[0]);
  sz2 = sq*sqrt(u[1]);
  sz3 = sq*sqrt(u[2]);

  z[0] = -0.25*b + sz1 + sz2 + sz3;
  z[1] = -0.25*b + sz1 - sz2 - sz3;
  z[2] = -0.25*b - sz1 + sz2 - sz3;
  z[3] = -0.25*b - sz1 - sz2 + sz3;
  SortArray(z,4);

  if (debug_print){
    printLog ("Quartic roots = %f  %f  %f  %f; q = %8.3e\n",z[0],z[1],z[2],z[3],q);
    CheckSolutions(b,c,d,e,z);
  }

  return 0;
}

/* ********************************************************************* */
void CheckSolutions(double b, double c, double d, double e, double *z)
/*
 *
 ********************************************************************** */
{
  int n;
  double x, f, fp, fm, fmax;

/* -- scale to min and max solution */

  x  = 1.0;  fp = e + x*(d + x*(c + x*(b + x)));
  x  = -1.0; fm = e + x*(d + x*(c + x*(b + x)));
  for (n = 0; n < 4; n++){
    x = z[n];
    f =  e + x*(d + x*(c + x*(b + x)));
    printLog ("> CheckSolutions(): f(z[%d]) = %8.3e\n",n,f);
  }
   
}

/* ********************************************************************* */
void QuarticPrintCoeffs(double b, double c, double d, double e)
/*
 *
 ********************************************************************** */
{
  printLog ("  f(x) = %18.12e + x*(%18.12e + x*(%18.12e ",
         e, d, c);
      printLog ("  + x*(%18.12e + x*%18.12e)))\n", b, 1.0);

  printLog ("  b = %18.14e;\n",b);
  printLog ("  c = %18.14e;\n",c);
  printLog ("  d = %18.14e;\n",d);
  printLog ("  e = %18.14e;\n",e);

  double a2 = -2.0*c;
  double a1 = c*c + b*d - 4.0*e;
  double a0 = -(b*c*d - b*b*e - d*d);

  printLog ("  Resolvent cubic:\n");
  printLog ("  g(x) = %18.12e + x*(%18.12e + x*(%18.12e + x))\n", a0, a1, a2);
  double Y    = 0.25*b*b + 2.0*d/(b+1.e-40) - c;
  printLog ("  b^2/4 = %f\n",0.25*b*b);
  printLog ("  Y     = %8.3e\n",Y);  
}


/* ********************************************************************* */
void PrintSolution (double *z)
/*
 *
 ********************************************************************** */
{
  printLog ("z = %f  %f  %f  %f\n",z[0],z[1],z[2],z[3]);
}



int CubicSolve2 (double b, double c, double d, double z[])
{
  double Q, R, R_Q, B, A, sQ, th, tn, sn, cn;

  Q  = (b*b - 3.0*c)/9.0;                       /* Eq. [5.6.10] */
  R  = (2.0*b*b*b - 9.0*b*c + 27.0*d)/54.0;  /* Eq. [5.6.10] */
  R_Q = R/Q;

if (debug_print) printLog ("R = %8.3e, Q = %8.3e, R^2/Q^3 - 1.0 = %8.3e\n",R,Q,R_Q*R_Q/Q-1.0);
  if (R_Q*R_Q <= Q){   /*** First case: 3 real roots ***/
if (debug_print) printLog ("> CubicSolve2(): 3 real roots\n");    
    sQ = sqrt(Q);
    th = acos(R_Q/sQ);    /* Eq. [5.6.11]  */
    tn = tan(th/3.0);
    cn = 1.0/sqrt(1.0 + tn*tn);
    sn = cn*tn;

  /* The three roots are given by Eq. (5.6.12) of Numerical Recipe
     (we switch x3 <-> x2 so that one can verify that x1 < x2 < x3 always): 

     x1 = -2.0*sQ*cos(th/3.0) - aa/3.0;
     x2 = -2.0*sQ*cos((th - 2.0*CONST_PI)/3.0) - aa/3.0;   
     x3 = -2.0*sQ*cos((th + 2.0*CONST_PI)/3.0) - aa/3.0;

     To avoid loss of precsion we use
       
       cos((th + 2pi/3) = cos(th/3)*cos(2pi/3) - sin(th/3)*sin(2pi/3) 
       cos((th - 2pi/3) = cos(th/3)*cos(2pi/3) + sin(th/3)*sin(2pi/3) 

     Using the fact that the roots must lie in [0,1] and that the spline
     is monotonically increasing, this is how we pick up the solution:   */ 


     z[0] = -2.0*sQ*cos(th/3.0) - 1.0/3.0;
     z[1] = -2.0*sQ*cos((th - 2.0*CONST_PI)/3.0) - b/3.0;   
     z[2] = -2.0*sQ*cos((th + 2.0*CONST_PI)/3.0) - b/3.0;
      
  }else{    /*** Second case: one root only ***/
    printLog( "! CubicSolve2(): 1 root only, arg = %8.3e !!\n", R_Q*R_Q/Q);

    A = Q/R;
    A = fabs(R)*(1.0 + sqrt(1.0 - A*A*Q));
    A = -DSIGN(R)*pow(A, 1.0/3.0);    /* Eq. [5.6.15] */
    if (A == 0) B = 0.0;
    else        B = Q/A;
    z[0] = (A + B) - b/3.0;            /* Eq. [5.6.17] */
    z[1] = z[2] = 0.0;
  }


  if(debug_print) {
    int l;
    double x,f;
    printLog ("===========================================================\n");
    printLog ("> Resolvent cubic:\n");
    printLog ("  g(x) = %18.12e + x*(%18.12e + x*(%18.12e + x))\n", d, c, b);
    printLog ("> Cubic roots = %8.3e  %8.3e  %8.3e\n",z[0],z[1],z[2]);
    for (l = 0; l < 3; l++){  // check accuracy of solution
     
      x = z[l];
      f = d + x*(c + x*(b + x));
      printLog ("  verify: g(x[%d]) = %8.3e\n",l,f);   
    }

    printLog ("===========================================================\n");
  }

  return 0;
}
