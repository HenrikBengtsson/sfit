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
 * cfit.c - implementation of simplex-fitting algorithm
 *
 */
#include "cfit.h"

void
cfit_dump ( struct cfit_context *C )
{
  if( C->outX )
    A_out ( C->n, C->k, C->X, C->outX ); /* BUG: if proj-xform is used, the
                                         coefficients is in the projection
                                         space. */
  
  if( C->outM )
    {
    if( C->gamma <= 0 )
      A_out ( C->k, C->m, C->M, C->outM );
    else /* output inverse-projective transformed simplex. */
         /* Note that this is NOT the simplex used in iteration */
         /* WARNING! Mnext is clobbered! */
      {
      int i;
      for ( i = 1; i <= C->k; i++ )
        gvFu_w ( C->m, C->gamma, C->M[i], C->u, C->Mnext[i] );
      if( C->o )
        Ap1vT ( C->k, C->m, C->Mnext, C->o );
      A_out ( C->k, C->m, C->Mnext, C->outM );
      }
    }
  
  if( C->gamma > 0 && C->outPM ) /* the simplex used in iteration */
    A_out ( C->k, C->m, C->M, C->outPM );
}

/* The simplex-fitting algorithm
 *
 * The weird reordering of formulas and loops is meant to optimize
 * memory access (minimize cache miss). Large matrices are accesses
 * in a single-pass and sequentially (no large strides) if possible.
 * The assumption is n > m >> k. Multiple cache lines is assumed.
 * 
 * Presently, there is no benchmark to support that this scheme is better
 * than a more straightforward implementation.
 * 
 */

void
cfit ( struct cfit_context *C )
{
  int h, i, j, iter, ia, rank;
  int n = C->n, m = C->m, k = C->k;
  R *x, *mu, *y, *th, *sth, *w, *rho;
  R alpha, frobratio; /* ratio of frobenius norm of DM and M */
#ifdef WITH_ROBUST
  R beta;
#endif
  
  for ( iter = 0, ia = 1; ia <= C->nalpha && iter < C->maxiter; ia++ )
    {
    alpha = C->alpha[ia];
#ifdef WITH_ROBUST
    beta = - C->robustcoeff * alpha;
#endif
    if( C->dump_mode == EVERY_ALPHA )
      cfit_dump ( C );
    
    for( frobratio = R_MAX; frobratio > C->frobratio && iter < C->maxiter; )
      {
      if ( C->dump_mode == EVERY_REX )
        cfit_dump ( C );

      /* R-step */
      Am1vT_B ( k - 1, m, C->M, C->M[k], C->Mi ); 
      Am1vT_B ( n, m, C->Y, C->M[k], C->Yi );
      rank = qrpsolve ( m, k - 1, n, C->Mi, C->Yi, C->X,
          C->QRtolerance, k, 0, C->normq, C->permut );
      
      /* make X[j] affine, compute tau_ij weight_j -> W_ij, and sum_j W_ij */
      c_v ( k, 0, C->sumWj );
      for( x = C->X[1], w = C->W[1], j = 1; j <= n; j++, x += k, w += k )
        {
        R phi, sumphi = 0, sumx = 0;
        for ( i = 1; i < k; i++ )
          {
          phi = x[i] < 0 ? 0 : ( x[i] > 1 ? 1: x[i] );
          w[i] = phi * C->weight[j];
          sumphi += phi;
          sumx += x[i];
          }
        x[i] = 1 - sumx;
        phi = x[i] < 0 ? 0 : ( x[i] > 1 ? 1: x[i] );
        w[i] = phi * C->weight[j];
        sumphi += phi;
        for ( i = 1; i <= k; i++ )
          {
          w[i] /= sumphi;       /* sumphi == 0 should be impossible */
          C->sumWj[i] += w[i];
          }
        }
      
      /* compute residual */
      for ( y = C->Y[1], rho = C->Rho[1], x = C->X[1], j = 1;
            j <= n; j++, y += m, rho += m, x += k )
        {
        u_v ( m, y, rho );
        for ( mu = C->M[1], i = 1; i <= k; i++, mu += m )
          cupv ( m, -x[i], mu, rho );
        }
      
      /*
       * X-step: delta mu_i due to residuals
       */
      /* accumulate nominator (denominator already computed above) */
      c_v ( k * m, 0, C->Mnext[1] );
      for ( w = C->W[1], j = 1, rho = C->Rho[1]; j <= n; j++, w += k, rho += m )
        for ( mu = C->Mnext[1], i = 1; i <= k; i++, mu += m )
          if( C->normq[i] == 0 )
            continue; /* do not update if basis is unused */
          else
            cupv ( m, w[i], rho, mu );
      /* normalize */
      for ( mu = C->Mnext[1], i = 1; i <= k; i++, mu += m )
        if( C->normq[i] > 0 && C->sumWj[i] > 0)
          cv ( m, 1/C->sumWj[i], mu );
      
      /*
       * X-step: delta mu_i due to theta_{hi}
       */
      c_v ( k * k, 0, C->Theta[1] );
      c_v ( k * k, 0, C->sumWPsi[1] );
      
      /* accumulate the numerator and denominator */
      for ( w = C->W[1], x = C->X[1],
          j = 1; j <= n; j++, w += k, x += k )
        for ( sth = C->sumWPsi[1], th = C->Theta[1],
            i = 1; i <= k; i++, sth += k, th += k )
          for ( h = 1; h <= k; h++ )
            {
            R f;
            if ( h == i ) continue;
#ifdef WITH_ROBUST
            if( x[h] > 0 )
              { f = w[i] * alpha; th[h] += f * x[h]; sth[h] += f; }
            else if( x[h] < beta )
              { th[h] += w[i] * beta; sth[h] += beta / x[h]; }
            else
              { f = w[i] * (1-alpha); th[h] += x[h]; sth[h] += f; }
#else /* expectile */            
            f = w[i] * ( x[h] < 0 ? 1 - alpha : alpha );
            sth[h] += f;
            th[h] += f * x[h];
#endif
            }
      
      /* normalize */
      /*
      for ( th = C->Theta[1], sth = C->sumWPsi[1], i = 1; i <= k * k; i++ )
        if( sth[i] != 0 )
          th[i] /= sth[i];
      */
      for ( th = C->Theta[1], sth = C->sumWPsi[1], i = 1;
          i <= k; i++, th += k, sth += k )
        {
        R f = 0;
        for ( h = 1; h <= k; h++ )
          if( h == i ) continue;
          else
            {
            if( sth[h] > 0 ) th[h]/= sth[h];
            f += th[h];
            }
        th[i] = 1 - f;
        }

      /* transform */
      for ( th = C->Theta[1], mu = C->Mnext[1], i = 1; i <= k; i++, th += k, mu += m )
        for ( x = C->M[1], h = 1; h <= k; h++, x += m)
          cupv ( m, th[h], x, mu );

      
      /* apply constraints if any, M and Mnext is modified */
      if( C->constraints )
        C->cf_function ( k, m, C->M, C->Mnext, C->cf_parameters);

      frobratio = vmuN22 ( k * m, C->Mnext[1], C->M[1] )
                  / vTv ( k * m, C->M[1] );
      u_v ( k * m, C->Mnext[1], C->M[1] ); /* make new vertices current */
      
      iter++;
      if( C->verbose)
        {
        R norm_xmin = 0, norm_xplus = 0, normRho;
        for ( x = C->X[1], i = 1; i <= n * k; i++ )
          if( x[i] < 0 ) norm_xmin += x[i]*x[i];
          else norm_xplus += x[i]*x[i];
        normRho = vTv ( n * m, C->Rho[1] );

        fprintf(stderr,"%3d  %8.6"Rf"  %10"Re"  %10.8"Rf"  %10"Re" %c\n",
            iter, alpha, frobratio,
            norm_xmin/(norm_xmin+norm_xplus), normRho,
            rank < k-1 ? '*':' ');
        }
      } /* REX */
    } /* alpha */
}

/*
 * Find an initial simplex based on incomplete QR factorization
 * of the data.
 */
void
initial_simplex ( struct cfit_context *C )
{
  int i, j, jmax; R *y, *yi, *mu, distmax;
  int n = C->n, m = C->m, k = C->k;
  int *P = new_P ( n );
  
  mu = C->M[k];
  c_v ( m, 0, mu );
  for ( y = C->Y[1], j = 1; j <= n; j++, y += m )
    cupv ( m, C->weight[j], y, mu );
  for ( distmax = -R_MAX, jmax =0, y = C->Y[1], j = 1; j <= n; j++, y+= m )
    {
    R dist = C->weight[j] * vmuN22 ( m, y, mu );
    if( dist > distmax )
      { jmax = j; distmax = dist; }
    }
  u_v ( m, C->Y[jmax], mu );
  for ( y = C->Y[1], yi = C->Yi[1], j = 1; j <= n; j++, y+=m, yi+=m )
    for ( i = 1; i <= m; i++ )
      yi[i] = C->weight[j] * ( y[i] - mu[i] );

  qrpsolve ( m, n, n, C->Yi, NULL, C->Y,
          C->QRtolerance, 0, k - 1, NULL, P );

  for ( i = 1; i < k; i++ )
    u_v ( m, C->Y[ P[i] ], C->M[i] );

  i_del ( P );
  return;
} 
  
/* projecting to the wall of positive orthant is fairly trivial */
void cf_nn ( int k, int m, R **M, R **Mnext, void *param)
{
  int i;
  for ( i = 1; i <= k*m; i++ )
    if( Mnext[i] < 0 ) Mnext[i] = 0; 
}

/* we assume the parameters are pointers to 2 m-vectors
 * corresponding to the lower and upper bound, respectively
 */
void cf_box ( int k, int m, R **M, R **Mnext, void *param)
{
  int i, j;
  R *l, *u;

  if( !param ) return;
  l = (R*)param; u = ((R*)param) + m;
  for ( i = 1; i <= k; i++ )
    for ( j = 1; j <= m; j++ )
      {
      if( Mnext[i][j] < l[j] ) Mnext[i][j] = l[j];
      if( Mnext[i][j] > u[j] ) Mnext[i][j] = u[j];
      }
}

/* Probability-like constraints
 * 
 * The algorithm here is based on Luenberger's gradient projection
 * (section 11.4, p332-333). The gradient is the distance between
 * the current vertex and the next one.
 *
 * The constraints here is very special:
 *  I mu >= 0 and 1^T mu = 1
 * The expression -(A A^T)^-1 A grad f^T  and
 * P grad =  { I - A^T(A A^T) A } grad f^T
 * can be computed easily.
 *
 * Note that for simplicity grad f has nothing to do with
 * the simplex fitting algorithm, we just assume the closest
 * Euclidean distance to the probability simplex (the feasible set)
 * is what we want to find. In another words, it's a quadratic
 * program that minimizes || mu - mu0 ||_2, i.e. the Hessian is I.
 */
#define TINY sqrt(R_EPSILON)
#define CF_PROB_MAXITER 30
void cf_prob ( int k, int m, R **M, R **Mnext, void *param)
{
  int i, j, iter, *active = new_i ( m );
  R *mu, *mu0, *d = new_v ( m ), *mutry = new_v ( m );
  for ( mu = Mnext[1], mu0 = M[1],
      i = 1; i <= k; i++, mu += m, mu0 += m )
    {
    for ( iter = 1; iter <= CF_PROB_MAXITER; iter++)
      {
      R sumgi, normd, dmax, dmin; int jmin, jmax;
      int n; /* number of inactive */

      /* find active set */
      for ( sumgi = 0, n = 0,  j = 1; j <= m; j++ )
        if ( fabs(mu0[j]) < TINY )
          { active[j] = 1; }
        else
          { active[j] = 0; n++; sumgi += mu[j] - mu0[j]; }
      sumgi /= n;
      
      FIND_DIRECTION:
      for ( normd = 0, dmax = -R_MAX, jmax = 0, j = 1; j <= m; j++ )
        {
        d[j] = active[j] ? 0 : mu[j] - mu0[j] - sumgi;
        normd += d[j] * d[j];
        if( d[j] > dmax ) { dmax = d[j]; jmax = j; }
        }
      if ( normd < TINY ) /* d parallel to the active subspace */
        {
        if( dmax <= TINY )
          break;            /* KKT condition applies, we're done */
        
        /* otherwise, delete constraint with the largest lambda */
        active[jmax] = 0;
        n--;
        if( n == 0 )
          sumgi = 0;
        else 
          sumgi = (sumgi * (n+1) - (mu[jmax]-mu0[jmax]))/n;
        goto FIND_DIRECTION;
        }
      else
        {
        R alpha;
        vpu_w ( m, d, mu0, mutry );
        for( jmin = 0, dmin = R_MAX, j = 1; j <= m; j++ )
          if( mutry[j] < dmin ){ dmin = mutry[j]; jmin = j; }
        if( dmin < -TINY )
          {
          alpha = - mu0[jmin] / d[jmin];
          cupv ( m, alpha, d, mu0 );
          }
        else
          u_v ( m, mutry, mu0 );
        }
      }
      u_v ( m, mu0, mu );
    }
  i_del ( active );
  v_del ( d );
  v_del ( mutry );
}
