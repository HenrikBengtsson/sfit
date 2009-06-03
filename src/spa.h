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
 * spa.h - an ad hoc (but no frill) signal-processing library
 *
 * - define vector and matrix data structure conventions
 * - mostly inline'd functions (should work with gcc or C99 compilers)
 * - some basic operators
 * - not meant to be fast (you should look at ATLAS/BLAS for that), but
 *   only to be as fast as naive codes can be. The main principle is array
 *   access should be as contiguous as possible, given row-major order.
 * 
 */
#ifndef _SPA_H_
#define _SPA_H_

#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Use the type R for floating point throughout. The macro Re Rf and Rg
 * are used in the format strings of scanf and printf, eg. "%9.7"Re"\n".
 *
 */
#ifdef R_IS_FLOAT  /* single precision */
typedef float R;
#define R_EPSILON FLT_EPSILON
#define R_MAX FLT_MAX
#define R_MIN FLT_MIN
#define Re "e"
#define Rf "f"
#define Rg "g"

#else /* double precision */
typedef double R;
#define R_EPSILON DBL_EPSILON
#define R_MAX DBL_MAX
#define R_MIN DBL_MIN
#define Re "le"
#define Rf "lf"
#define Rg "lg"

#endif /* choice of precision */

#define IN_VEC_CHUNK 64    /* reallocation chunks for autoallocated input */
#define IN_MAT_CHUNK 256


/* Vectors and matrices are unit-indexed. Matrices are row-major and
 * allocated as a contiguous block. An array of row vector is provided
 * for random access. Note that they point to possibly non-existent zero-th
 * element.
 *
 * We adopt a `row-centric' view here, unlike standard linear algebra notation.
 * This is because in practice, the interesting vectors are more naturally
 * packaged (e.g. in data files or display) as row vectors. Consequently,
 * operations like linear transformation is done by right-multiplication.
 * 
 * Function names are terse. We use a mnemonic convention that hopefully
 * describe what a function does succintly. v stands for vectors and A
 * for matrix. A function such as vA_u means right-multiply v by A and put
 * the result in u.
 */


/** CONSTRUCTOR AND DESTRUCTOR **/

static inline R*
new_v ( int m )                  /* return a new m-vector */
{
  R* v = (R*) malloc ( sizeof(R) * m );
  return v ?  v - 1 : NULL;
}
  
static inline void
v_del ( R* x )                 /* free a vector */
{
  free( & x[1] );
}

static inline R**
new_A ( int n, int m )          /* return an nxm matrix */
{
  int i;
  R** A = (R**) malloc ( sizeof(R*) * n );
  if( A ) A-- ; else return NULL;
  A[1] = (R*) malloc ( sizeof(R) * n * m );
  if( A[1] ) A[1]-- ; else { free( &A[1] ); return NULL; }
  for ( i = 2; i <= n; i++ )
    A[i] = A[i-1] + m;
  return A;
}

static inline void
A_del ( R** A )                /* free a matrix **/
{
  free( &A[1][1] ); free( &A[1] );
}

static inline int*
new_i ( int m )                /* integer vector */
{
  int *i = (int*) malloc ( sizeof(int) * m );
  return i ? i - 1 : NULL;
}

static inline int* 
new_P ( int m )                /* permutation, int vector filled with 1,...,m */
{
  int i, *P = (int*) malloc ( sizeof(int) * m );
  if( P ) P-- ; else return NULL;
  for ( i = 1; i <= m; i++ ) P[i] = i;
  return P;
}

static inline void
i_del (int *i)    /* use this to free permutation as well */
{
  free( &i[1] );
}

/** INPUT AND OUTPUT **/

static inline int
in_v ( FILE *in, int m, R* v ) /* scan m numbers from in */
{
  int i;
  for ( i = 1; i <= m; i++ )
    if( 1 != fscanf( in, " %"Rf, & v[i] ) ) break;
  return i - 1;                         /* the number successfully read */
}

extern R*
in_newv ( FILE *in, int *m );
  /* Read [ \t\v\r\f] separated numbers until the end-of-line.
     Return newly allocated vectors, m is the number of elements found.
     realloc() may be used repeatedly, the macro IN_VEC_CHUNK specify how
     grainy the chunks.
   */

static inline void
v_out ( int m, R* v, FILE *out ) /* print m numbers to out */
{
  int i;
  for ( i = 1; i < m; i++ )
    fprintf( out, "%"Rg"\t", v[i] );
  fprintf( out, "%"Rg"\n", v[i] );
}

static inline void
v_fout ( int m, R* v, char *format, FILE *out ) /* print with custom format */
{
  int i;
  for ( i = 1; i <= m; i++ )
    fprintf( out, format, v[i] );           /* users supply the separator ! */
  fputc('\n', out);
}

static inline void
i_out ( int m, int *v, FILE *out )
{
  int i;
  for ( i = 1; i < m; i++ )
    fprintf( out, "%d\t", v[i] );
  fprintf( out, "%d\n", v[i] );
}

static inline int
in_A ( FILE *in, int n, int m, R** A )   /* read matrix data */
{
  int i, j; R *a; 
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      if( 1 != fscanf(in, " %"Rf, & a[j] ) ) return (j - 1) + (i-1)*m;
  return n * m;
}

extern R**
in_newA ( FILE *in, int *n, int *m ); /* read and auto-allocate matrix */

static inline void
A_out ( int n, int m, R **A, FILE *out )
{
  R* a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    {
    for ( j = 1; j < m; j++ )
      fprintf( out, "%"Rg"\t", a[j] );
    fprintf( out, "%"Rg"\n",a[j]);
    }
}

static inline void
A_fout ( int n, int m, R **A, char* format, FILE *out )
{
  R* a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    {
    for ( j = 1; j <= m; j++ )
      fprintf( out, format, a[j] );
    fputc('\n', out );
    }
}

/*
 *    FUNCTIONALS
 */
static inline R
uTv ( int m, R *u, R *v )   /* u^T v = (u,v) = u.v ,  dot product */
{
  R f; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    f += u[i] * v[i];
  return f;
}

static inline R
vTv ( int m, R *v )  /* v^T x = || v ||_2^2 */
{
  R f; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    f += v[i] * v[i];
  return f;
}

static inline R
vT1 ( int m, R *v )          /* v^T 1 = sum v_i */
{
  R f; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    f += v[i];
  return f;
}

static inline R             /* 1-norm,  sum |v_i| */
vN1 ( int m, R *v )
{
  R f; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    f += v[i] < 0 ? -v[i] : v[i];
  return f;
}

static inline R
vN2 ( int m, R *v )         /* 2-norm, sqrt ( sum v_i^2 ) */
{
  R f; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    f += v[i] * v[i];
  return sqrt(f);
}

static inline R
vmuN1 ( int m, R *v, R *u )  /* || v - u ||_1 */
{
  R f, c; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    {
    c = v[i] - u[i];
    f += c < 0 ? -c : c;
    }
  return f;
}

static inline R
vmuN22 ( int m, R *v, R *u ) /* || v - u ||_2^2   */
{
  R f, c; int i;
  for ( f = 0, i = 1; i <= m; i++ )
    {
    c = u[i] - v[i];
    f += c * c;
    }
  return f;
}

static inline int
imax ( int m, R *v )          /* find the first maximum */
{
  R max = -R_MAX; int i, imax = 0;
  for ( i = 1; i <= m; i++ )
    if( v[i] > max ) { imax = i; max = v[i]; }
  return imax;
}

static inline int
imin ( int m, R *v )          /* find the first maximum */
{
  R min = R_MAX; int i, imin = 0;
  for ( i = 1; i <= m; i++ )
    if( v[i] > min ) { imin = i; min = v[i]; }
  return imin;
}

/*  OPERATION ON THE ELEMENTS
 *  
 *  Because of the continuous block in matrices, matrix element-wise operation
 *  can be done by considering A[1] a vector of length n*m. For example,
 *  A - B -> C is done using  vmu_w ( n * m, A[1], B[1], C[1] ).
 */
static inline void
u_v ( int m, R *u, R* v )     /* copy u to v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] = u[i];
}

static inline void
uv_vu ( int m, R *u, R* v )   /* swap u and v */
{
  int i; R c;
  for ( i = 1; i <= m; i++ )
    { c = v[i]; v[i] = u[i]; u[i] = c; }
}

static inline void
c_v ( int m, R c, R *v )  /* v = ( c, c, c, ..., c ) */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] = c;
}

static inline void
cmore_v ( int m, R c, R *v )   /* if v[i] > c, set to c */
{
  int i;
  for ( i = 1; i <= m; i++ )
    if( v[i] > c ) v[i] = c;
}

static inline void
cless_v ( int m, R c, R *v )    /* if v[i] < c, set to c */
{
  int i;
  for ( i = 1; i <= m; i++ )
    if( v[i] < c ) v[i] = c;
}

static inline void
dirac (int m, int j, R *v )     /* delta vector */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] = 0;
  v[j] = 1;
}
  
/* function applications (maps) */
static inline void
fv ( int m, double (*f)(double), R *v )    /* f(v[i]) -> v[i] */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] = f(v[i]);
}

static inline void
fv_u ( int m, double (*f)(double), R *v, R *u )    /* f(v[i]) -> u[i] */
{
  int i;
  for ( i = 1; i <= m; i++ )
    u[i] = f(v[i]);
}

static inline void
frange_v ( int m, R v1, R dv, double (*f)(double), R *v )
{                                         /* f( v1 + (i-1) dv) -> v[i] */
  int i;
  if( !f ) /* f is identity */
    for ( i = 1; i <= m; i++ )
      v[i] = v1 + (i-1) * dv;
  else
    for ( i =1; i <= m; i++ )
      v[i] = f( v1 + (i-1) * dv );
}

static inline void
rint_v ( int m, R *v )  /* round to integer */
{
  int i;
  for ( i = 0; i <= m; i++ )
    v[i] = rint(v[i]);
}

static inline void
round_v ( int m, R resolution, R *v ) /* round to the nearest multiple of */
{                                     /* `resolution' */
  int i;
  for ( i = 0; i <= m; i++ )
    v[i] = resolution * rint(v[i]/resolution);
}

/* scalar-vector ops */

static inline void
cv ( int m, R c, R *v )       /* c v -> v,    multiply v by c, in place */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] *= c;
}

static inline void
cv_u ( int m, R c, R *v, R *u ) /* c v -> u */
{
  int i;
  for ( i = 1; i <= m; i++ )
    u[i] = c * v[i];
}

static inline void
opv ( int m, R o, R *v )     /* o + v -> v,    offset by a constant */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] += o;
}

static inline void
opv_u ( int m, R o, R *v, R *u )   /* o + v -> u */
{
  int i;
  for ( i = 1; i <= m; i++ )
    u[i] = o + v[i];
}

static inline void
opcv ( int m, R o, R c, R *v )     /* o + c v -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] = o + c * v[i];
}

static inline void
opcv_u ( int m, R o, R c, R *v, R *u )     /* o + c v -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    u[i] = o + c * v[i];
}

/* vector-vector ops */
static inline void
vpu_v ( int m, R *v, R *u )   /* u + v -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] += u[i];
}

static inline void
vpu_w ( int m, R *v, R *u, R *w )  /* u + v -> w */
{
  int i;
  for ( i = 1; i <= m; i++ )
    w[i] = u[i] + v[i];
}

static inline void
vtu_v ( int m, R *v, R *u )   /* u v -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] *= u[i];
}

static inline void
vtu_w ( int m, R *v, R *u, R *w )  /* u v -> w */
{
  int i;
  for ( i = 1; i <= m; i++ )
    w[i] = u[i] * v[i];
}

static inline void
vmu_v ( int m, R *v, R *u )   /* v - u -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] -= u[i];
}

static inline void
vmu_w ( int m, R *v, R *u, R *w )  /* v - u -> w */
{
  int i;
  for ( i = 1; i <= m; i++ )
    w[i] = v[i] - u[i];
}

static inline void
vdu_v ( int m, R *v, R *u )   /* v / u -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] /= u[i];
}

static inline void
vdu_w ( int m, R *v, R *u, R *w )  /* v / u -> w */
{
  int i;
  for ( i = 1; i <= m; i++ )
    w[i] = v[i] / u[i];
}

static inline void
cupv ( int m, R c, R *u, R *v )    /* c u + v -> v */
{
  int i;
  for ( i = 1; i <= m; i++ )
    v[i] += c * u[i];
}

static inline void
cupv_w ( int m, R c, R *u, R *v, R *w )    /* c u + v -> w */
{
  int i;
  for ( i = 1; i <= m; i++ )
    w[i] = c * u[i] + v[i];
}

static inline R
vF1 ( int m, R *v )        /* v/(1^T v) -> v,   return 1^T v */
{
  R F; int i;
  for ( F = 0, i = 1; i <= m; i++ )
    F += v[i];
  for ( i = 1; i <= m; i++ )
    v[i] /= F;
  return F;
}

static inline R
vF1_u ( int m, R *v, R *u )        /* v/(1^T v) -> u,   return 1^T v */
{
  R F; int i;
  for ( F = 0, i = 1; i <= m; i++ )
    F += v[i];
  for ( i = 1; i <= m; i++ )
    u[i] = v[i] / F;
  return F;
}

/* projective transform & inverse */
static inline R
vFug ( int m, R *v, R *u, R g )      /*   v/(u^T v + g ) -> v */
{
  R F; int i;
  for ( F = 0, i = 1; i <= m; i++ )
    F += u[i] * v[i];
  F += g;
  for ( i = 1; i <= m; i++ )
    v[i] /= F;
  return F;
}

static inline R
vFug_w ( int m, R *v, R *u, R g, R *w )      /*   v/(u^T v + g ) -> w */
{
  R F; int i;
  for ( F = 0, i = 1; i <= m; i++ )
    F += u[i] * v[i];
  F += g;
  for ( i = 1; i <= m; i++ )
    w[i] = v[i] / F;
  return F;
}

static inline R
gvFu ( int m, R g, R *v, R *u )     /* g v / ( 1 - u^T v ) -> v */
{
  R F; int i;
  for ( F = 0, i = 1; i <= m; i++ )
    F += u[i] * v[i];
  F = (1 - F)/g;
  for ( i = 1; i <= m; i++ )
    v[i] /= F;
  return F;
}

static inline R
gvFu_w ( int m, R g, R *v, R *u, R *w )     /* g v / ( 1 - u^T v ) -> w */
{
  R F; int i;
  for ( F = 0, i = 1; i <= m; i++ )
    F += u[i] * v[i];
  F = (1 - F)/g;
  for ( i = 1; i <= m; i++ )
    w[i] = v[i] / F;
  return F;
}

/*
 * MATRIX OPERATIONS
 *
 */

static inline void  /* transpose (and invert) permutation matrix */
P_PT ( int m, int *P, int *PT )
{
  int i;
  for ( i = 1; i <= m; i++ )
    PT[P[i]] = i;
}

static inline void
AT ( int n, R **A )   /* in-place transpose, A has to be square */
{
  int i, j; R c, *a, *b;
  for ( a = A[1], i = 1; i < n; i++, a += n )
    for ( b = a + n + i, j = i+1; j <= n; j++, b += n )
      { c = *b; *b = a[j]; a[j] = c; }
}

static inline void
AT_B ( int n, int m, R **A, R **B )   /* A^T -> B */
{
  int i, j; R *a = A[1], *b = B[1];
  for ( i = 1; i <= n; i++, a += m )
    for ( b = &B[1][i], j = 1; j <= m; j++, b += n )
      *b = a[j];
}

static inline void
AB_C ( int ni, int nj, int nk, R **A, R **B, R **C ) /* matrix multiplication */
{
  R *a, *b, *c; int i, j, k;
  for ( a = A[1], c = C[1], i = 1; i <= ni; i++, a += nj, c += nk )
    {
    b = B[1];
    for ( k = 1; k <= nk; k++ )
      c[k] = a[1] * b[k];
    for ( j = 2; j <= nj; j++ )
      for ( b += nk, k = 1; k <= nk; k++ )
        c[k] += a[j] * b[k];
    }
}

static inline void
ATB_C ( int ni, int nj, int nk, R **A, R **B, R **C ) /* A is nj x ni */
{
  R *a, *b, *c; int i, j, k;

  for ( a = A[1], b = B[1], c = C[1], i = 1; i <= ni; i++, c += nk )
    for ( k = 1; k <= nk; k++ )
      c[k] = a[i] * b[k];
  for ( a += ni, b += nk, j = 2; j <= nj; j++, a += ni, b += nk )
    for ( c = C[1], i = 1; i <= ni; i++, c += nk )
      for ( k = 1; k <= nk; k++ )
        c[k] += a[i] * b[k];
}

static inline void
ABT_C ( int ni, int nj, int nk, R **A, R **B, R **C ) /* B is nk x nj */
{
  R *a, *b, *c; int i, j, k;
  for ( a = A[1], c = C[1], i = 1; i <= ni; i++, a += nj, c += nj )
    for ( b = B[1], k = 1; k <= nk; k++, b += nj )
      for ( c[k] = 0, j = 1; j <= nj; j++ )
        c[k] += a[j] * b[j];
}

/* TODO: ATBT_C, AB_CT, ATB_CT, ABT_CT, ATBT_CT  !
 * there is a regularity in this. the innermost index should be
 * the one seen as rows twice, the outermost as column twice [ conjecture ]
 *
 * AB_CT is the only one problematic.
 * */

/* The following fours are matrix vector multiplication */
static inline void
ATv_u ( int n, int m, R **A, R *v, R *u ) /* v^T A -> u^T */
{
  R *a; int i, j;
  for ( a = A[1], j = 1; j <= m; j++ )
    u[j] = v[1] * a[j];
  for ( a += m, i = 2; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      u[j] += v[i] * a[j];
}

static inline void
AT1_u ( int n, int m, R **A, R *u ) /* A^T 1,   sum rows */
{
  R *a; int i, j;
  a = A[1];
  for ( j = 1; j <= m; j++ )
    u[j] = a[j];
  for ( a += m, i = 2; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      u[j] += a[j];
}

static inline void
Av_u ( int n, int m,  R **A, R *v, R *u ) /* A v -> u */
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( u[i] = v[1]*a[1], j = 2; j <= m; j++ )
      u[i] += v[j] * a[j]; 
}

static inline void
A1_u ( int n, int m, R **A, R *u )  /* A 1,   sum column */
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a+=m )
    for ( u[i] = a[1], j = 1; j <= m; j++ )
      u[i] += a[j];
}

/* These are operations between all columns or row vectors and a vector.
 * e.g. for shifting the coordinate system
 */
static inline void
Apv1T ( int n, int m, R **A, R *v )  /* accumulate v to each column */
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] += v[i];
}

static inline void
Apv1T_B ( int n, int m, R **A, R *v, R **B )  /* add v to each column */
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] + v[i];
}

static inline void
Amv1T ( int n, int m, R **A, R *v )  /* accumulate -v to each column */
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] -= v[i];
}

static inline void
Amv1T_B ( int n, int m, R **A, R *v, R **B )  /* subtract v to each column */
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] - v[i];
}

static inline void
Ap1vT ( int n, int m, R **A, R *v )  /* as above but on rows */
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] += v[j];
}

static inline void
Ap1vT_B ( int n, int m, R **A, R *v, R **B )  
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] + v[j];
}

static inline void
Am1vT ( int n, int m, R **A, R *v )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] -= v[j];
}

static inline void
Am1vT_B ( int n, int m, R **A, R *v, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] - v[j];
}

/* Diagonal matrix (represented as vectors). This can be seen as an
 * operation similar to the above for multiplication. */
static inline void
DA ( int n, int m, R *D, R **A )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] *= D[i];
}

static inline void
DA_B ( int n, int m, R *D, R **A, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] * D[i];
}

static inline void
AD ( int n, int m, R **A, R *D )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] *= D[j];
}

static inline void
AD_B ( int n, int m, R **A, R *D, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] * D[j];
}

/* For the sake of completeness we'll do the same for division
 * (might be useful, eg in inverses)
 *
 * The mnemonic `rD' meant to be the reciprocals of D
 * */

static inline void
rDA ( int n, int m, R *D, R **A )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] /= D[i];
}

static inline void
rDA_B ( int n, int m, R *D, R **A, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] / D[i];
}

static inline void
ArD ( int n, int m, R **A, R *D )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] /= D[j];
}

static inline void
ArD_B ( int n, int m, R **A, R *D, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[1], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] / D[j];
}

/* outer-product update. note the similarity with Ap1vT and ilks */
static inline void
ApuvT ( int n, int m, R **A, R *u, R *v )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] += u[i] * v[j];
}

static inline void
ApuvT_B ( int n, int m, R **A, R *u, R *v, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[i], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] + u[i] * v[j];
}

static inline void
AmuvT ( int n, int m, R **A, R *u, R *v )
{
  R *a; int i, j;
  for ( a = A[1], i = 1; i <= n; i++, a += m )
    for ( j = 1; j <= m; j++ )
      a[j] -= u[i] * v[j];
}

static inline void
AmuvT_B ( int n, int m, R **A, R *u, R *v, R **B )
{
  R *a, *b; int i, j;
  for ( a = A[1], b = B[i], i = 1; i <= n; i++, a += m, b += m )
    for ( j = 1; j <= m; j++ )
      b[j] = a[j] - u[i] * v[j];
}

/* TODO:
 *
 * - The mess above can be simplified using macros (but not by inlining
 *   because inlined functions can't be passed on). A code generator
 *   that can process mnemonics would be cool.
 * - This header file is too big. It might need to be split into most
 *   common operations, and least common ones, so that not everything
 *   is included.
 */

/* least-square solve using Householder QR and pivoting */
extern int         /* return the number of effective basis */
qrpsolve (
    int m, int k, int n,  /* A = k x m, Y = n x m, X = n x k */
    R **A, R **Y, R **X,  /* Y and X can be NULL */
    R min_ratio,          /* minimum ||q[k]||_2/||q[i]||_2 */
    int kX,               /* if kX > k, the row stride of X */
    int kQ,               /* if 0 < kQ< k, force to use the first kQ only. */
    R *unormq, int *uP    /* user may optionally ask for the norm
                              and permutation */
    );

#endif /* _SPA_H_ */
