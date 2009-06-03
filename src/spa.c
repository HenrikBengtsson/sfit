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
 * spa.c - implementations of non-inlined functions in spa.h
 *
 */

#include "spa.h"
R* in_newv ( FILE *in, int *m ) /* auto-allocate, '\n' sep'd */
{
  int i, nchunk, chr;
  R *x = NULL, c;

  nchunk = 0;
  i = 0;
  while ( fscanf(in," %"Rf"%*[ \t\v\r\f]",&c) == 1 )
    {
    if( i == nchunk )
      x = (R*) realloc ( x, sizeof(R) * ( nchunk += IN_VEC_CHUNK ) );
    x[i] = c;
    i++;
    
    chr = fgetc(in);
    if( chr == '\n' || chr == EOF ) break;
    ungetc(chr,in);
    }
  *m = i;
  return x - 1;
}

R** in_newA ( FILE *in, int *pn, int *pm )
{
  int n, m, i, j, chr, nchunk;
  R *x, *p = NULL, c;
  R **A;
  
  n = 0; m = 0;
  i = 0;      /* counter of items */
  x = NULL;
  nchunk = 0;
  j = 0;
  
  /* skip comments */
  while ( (chr=fgetc(in)) == '#' )
    { fscanf(in,"%*[^\n]"); fgetc(in); }
  ungetc(chr,in);
  
  while ( fscanf(in," %"Rg"%*[ \t\r]",&c) == 1 )
    {
    if( i == nchunk )
      {
      x = (R*) realloc ( x, sizeof(R) * ( nchunk += IN_MAT_CHUNK ) );
      p = & x[i];
      }
    *p++ = c;
    i++; j++;
    
    chr = fgetc(in);
    if( chr == '\n' || chr == EOF )
      {
      if ( n == 0 )
        m = j;      /* we use the first one as the vec length */
      else if( j != m )
        {
        fprintf(stderr,"non uniform row lengths\n" );
        free(p);
        *pn = *pm = 0;
        return NULL;
        }
      n++;
      j = 0;
      while ( (chr=fgetc(in)) == '#' )
        { fscanf(in,"%*[^\n]"); fgetc(in); }
      }
    ungetc(chr,in);
    }
  
  /* put it in a matrix */
  A = (R**) malloc ( sizeof(R) * n ) - 1;
  A[1] = x - 1;
  for ( i = 2; i <= n; i++ )
    A[i] = A[i-1] + m;
      
  *pn = n;
  *pm = m;
  return A;
}

/*
 * solve: min || Y - A X ||_2 using householder QRP^T
 *
 * If q[i]/q[1] < min_ratio qN22[j] <- 0,  x[j][i] <- 0
 *
 * A is overwritten by R and v, first k of y[j] is overwritten by
 * the permuted x[j].
 * 
 * 2/(v^T v) and v[i][i] and the permutation is lost.
 *
 * TODO:
 * - some of the code reorderings to improve memory access
 *   performance are still based on superstition. need to benchmark.
 * - the rank-deficiency detection (the stopping rule) might be too simplistic.
 */
int         /* return the number of effective basis */
qrpsolve (
    int m, int k, int n,  /* A = k x m, Y = n x m, X = n x k */
    R **A, R **Y, R **X,  /* Y and X can be NULL */
    R min_ratio,          /* minimum ||q[k]||_2/||q[i]||_2 */
    int kX,               /* if kX > k, the row stride of X */
    int kQ,               /* if 0 < kQ< k, force to use the first kQ only. */
    R *unormq, int *uP    /* user may optionally ask for the norm
                              and permutation */
    )
{
  int h, i, ek, j, *P; R* a, *b, *normq;
  
  if ( kQ > k || kQ <= 0 ) kQ = k;
  if ( kX < k ) kX = k;
  
  /* setup */
  if( uP )
    { P = uP; for ( i =1; i <= k; i++ ) P[i] = i; }
  else P = new_P ( k );
  normq = new_v ( k );
  for ( a = A[1], i = 1; i <= k; i++, a += m )
    normq[i] = vTv ( m, a );

  /* reduction of the basis vectors and the right-hand sides */
  for ( a = A[1], i = 1; i <= kQ; i++, a += m )
    {
    R mu, sigma, v1, beta, aTb;
    int p = i; R max = normq[i];
    
    for ( h = i + 1; h <= k; h++ )
      if( normq[h] > max )
        { max = normq[h]; p = h; }

    if( max / normq[1] < min_ratio )
      {
      for ( h = i; h <= k; h++ )
        normq[h] = 0;
      break;
      }

    if( p != i )  /* swap pivot */
      {
      j = P[p]; P[p] = P[i]; P[i] = j;
      normq[p] = normq[i]; normq[i] = max;
      uv_vu ( m, A[p], A[i] );
      }
    
    /* householder vector for i */
    for ( sigma = 0, j = i + 1; j <= m; j++ )
      sigma += a[j] * a[j];
    if ( sigma == 0 )
      { v1 = beta = 0; }
    else
      {
      mu = sqrt ( a[i] * a[i] + sigma );
      v1 = a[i] < 0 ? a[i] - mu : - sigma /(a[i]+mu);
      beta = 2/(sigma + v1*v1);
      a[i] = mu;
      }

    /* transform basis vectors */
    for ( b = a + m, h = i + 1; h <= k; h++, b += m )
      {
      for ( aTb = v1 * b[i], j = i+1; j <= m; j++ )
        aTb += a[j] * b[j];
      aTb *= beta;
      for ( b[i] -= aTb * v1, j = i+1; j <= m; j++ )
        b[j] -= aTb * a[j];
      normq[h] -= b[i] * b[i];
      }
    /* transform the right-hand sides, if supplied. */
    if( Y )
      for ( b = Y[1], h = 1; h <= n; h++, b += m )
        {
        for ( aTb = v1 * b[i], j = i+1; j <= m; j++ )
          aTb += a[j] * b[j];
        aTb *= beta;
        for ( b[i] -= aTb * v1, j = i+1; j <= m; j++ )
          b[j] -= aTb * a[j];
        }
    }
  ek = i-1; /* effective k, the last usable basis vector */
  if( Y && X )
    {
    /* backsubstitution (in-place)
     *
     * R is copied to another array, to make the access
     * pattern more regular when backsubstituting
     */
    R *r = new_v ( ek * ( ek + 1 )/2 );
    for ( j = 1, i = ek; i > 0; i-- )
      {
      for ( h = ek; h >= i; h-- )
        r[j++] = A[h][i]; 
      }

    /* This loop accesses a, b, r, and P at the same time.
     * Hopefully this doesn't thrash the cache. If this happens, permuting
     * the solution should have its own loop. */
    for ( b = Y[1], a = X[1], j = 1; j <= n; j++, b += m, a += kX )
      {
      int l;
      for ( l = 1, i = ek; i > 0; i-- )
        {
        for ( h = ek; h > i; h-- )
          b[i] -= b[h] * r[l++];
        b[i] /= r[l++];
        a[P[i]] = b[i];
        }
      for ( i = ek + 1 ; i <= k; i++ )
        a[P[i]] = 0;
      }
    v_del ( r );
    }
    
  if( unormq ) /* copy permuted norm */
    for ( i = 1; i <= k; i++ )
      unormq[P[i]] = normq[i];
  v_del ( normq );
  if( !uP ) i_del ( P );
  
  return ek;
}
