/*
 * ---------------------------------------------------
 * FOR RESTRICTED USE ONLY. DO NOT REDISTRIBUTE.
 * ---------------------------------------------------
 * 
 * Pratyaksha J. Wirapati
 * email: wirapati@wehi.edu.au
 * 
 * @2002
 * Genetics and Bioinformatics Division
 * Walter and Eliza Hall Institute
 * PO Royal Melbourne Hospital
 * Parkville, VIC 3050
 * Australia
 */
/*
 * cfit.h - interface to cfit algorithm
 */
#ifndef _CFIT_H_
#define _CFIT_H_

#include "spa.h"

struct cfit_context /* a data set (a computation task) and program options */
{
  /* data matrix */
  int n;            /* number of observations */
  int m;            /* number of variables */
  int k;            /* order of simplex, maximum number of basis */

  R *weight;        /* n        observation weights */
  R **Y;            /* n x m    data matrix */
  R **M;            /* k x m    current vertices */
  R **X;            /* n x k    coeff matrix */
  R **Rho;          /* n x m    residual vectors */
  
  R *normq;         /* k        norm of basis. used to flag discarded ones */
  int *permut;      /* k        permutation of vertices according to the norm */
  
  /* work space */
  R **Mnext;        /* k x m   (tentative) new vertices */
  R **Yi, **Mi;     /* n x m, n x k     Y and M w.r.t. a vertex */
  R **W;            /* n x k   `closeness' of an observation to a vertex */
  R **Theta;        /* k x k,  quantile-like stat */
  R **sumWPsi;      /* k x k,  sum of weight in M-estimator, psi(x)/x */
  R *sumWj;         /* k,      sum of W_ij over j */
  
  /* iteration controls */
  R *alpha;     /* sequence of alpha */
  int nalpha;       /* number of alpha's in alpha */
  
  int maxiter;     /* maximum REX steps */
  R frobratio;     /* tolerance for no improvement, ||DM||/||M|| */
  R QRtolerance;   /* tolerance for discarding basis vector */

  /* output files */
  FILE *inY, *inM;  /* input data file and (optional) initial values */
  FILE *outM, *outX; /* output files, for each matrix */
  FILE *outnormq;   /* output the norm of q_i as well */
  enum {            /* Indicates when the matrices are dumped. */
    NO_DUMP = 0,        /* Applied only to those out? which are not NULL */
    FINAL_ONLY = 1,     
    EVERY_ALPHA = 2,    /* matrices are concatenated */
    EVERY_REX = 3
  } dump_mode;

  /* vertex constraints */
  enum {
    NO_CONSTRAINTS = 0,
    NONNEGATIVE = 1,
    BOX = 2,
    PROBABILITY = 3
  } constraints;
  
  /* callback function that modify Mnext to satisfy some constraints */
  void (*cf_function)( int k, int m, R **M, R **Mnext, void *p);
  void *cf_parameters; /* constraint function  parameters */
  
  /* projective transform */
  R gamma; /* if negative, no proj. transform */
  R *u;    /* direction vector, the default is 1^T */
  R *o;    /* origin of transform, default 0 */
  char *fn_PY; /* dump proj-transformed Y to this file */
  FILE *outPM; /* dump proj-transformed M to this file */
  int denominator_as_weight; /* set the point weights proportional to denom. */
  
  /* miscellaneous */
  int verbose; /* verbose dumping to stderr */

  R robustcoeff;
  R chopquantileL, chopquantileH; /* truncation of data point */
};

extern void initial_simplex ( struct cfit_context *C );
extern void cfit ( struct cfit_context *C );

/* constraint functions */

extern void cf_nn ( int k, int m, R **M, R **Mnext, void *p );
  /* nonnegative */

extern void cf_box ( int k, int m, R **M, R **Mnext, void *p );
  /* boxed */

extern void cf_prob ( int k, int m, R **M, R **Mnext, void *p );
  /* probability-like */

/* output some of the matrices */
extern void cfit_dump ( struct cfit_context *C );

#endif
