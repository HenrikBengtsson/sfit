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
 *
 * HISTORY:
 * 2006-05-16
 * o Updated to compile with gcc 4.0.2.
 *
 */
/*
 * cli.c - command line interface to cfit algorithm
 *
 * This requires GNU getopt_long() function [ usually included in gcc ]
 * 
 */

#include "cfit.h"
#include <string.h>

#include "getopt.h"

void usage(void)
{ 
  fprintf(stderr,"\
cfit [ options ] k [ data [ weight ]]\n\
\n\
Fit a simplex with k vertices to points in the file <data> (stdin\n\
is used if not supplied).\n\
\n\
<weight> is a file of point weights (if not supplied, w_j = 1/n).\n\
\n\
-h, --help      this help\n\
-v, --verbose   dump progress to stderr\n\
\n\
-i, --Minit <filename>    initial simplex\n\
-d, --dump <int>          dump mode (0 none, 1 final, 2 each alpha, 3 each REX)\n\
-M, --outM <filename>     where M is dumped [ stdout ]\n\
-X, --outX <filename>     where X is dumped\n\
\n\
-a, --alpha <num,num,...> comma- or space-separated sequence of alpha\n\
-T, --maxiter <int>       maximum number of REX steps\n\
-f, --frobratio <num>     ratio of ||DM||/||M|| (inner loop convergence)\n\
    --QRtolerance <num>   ratio of ||q_i||/||q_1|| (QR termination)\n\
\n\
-C, --constraints  <int>  0 none, 1 nonnegative, 2 box, 3 probability-like\n\
-x, --lubox <lb,ub>       lower and upper bound in box constraints\n\
                          (specifying this turns on --constraints 2). If this\n\
                          is not specified, the dimension-specific lower and\n\
                          upper bounds are read from a file `lubox' containing\n\
                          2 x m  matrix of lower- and upper-bound vectors.\n\
\n\
-o, --origin <filename>   origin vector, offset data by this vector [default 0]\n\
-g, --gamma <num>        `eye-height' in projective transform [default 0]\n\
-u, --viewdir <filename> `viewing direction' vector [ meanY/(meanY^T meanY) ]\n\
    --u_1T                u = (1,1,1,... )\n\
-D, --denom_weight        use the denominators as observation weights\n\
    --outPY <filename>    save the transformed Y into this file\n\
    --outPM <filename>    dump projective-transformed M (according\n\
                          to the dump mode option).\n\
\n\
-q, --chopless <num>      assign zero weight to points if u^T y is less\n\
-Q, --chopmore <num>      or more than the specified percentiles.\n\
\n\
-R, --robust <num>        robust coefficient (if compiled with -DWITH_ROBUST)\n\
");
  exit(1);
}

/* callback for qsort */
int cmpR_asc ( const void *a, const void *b )
{
  return *(R*)a < *(R*)b ? -1 : 1;
}

int main( int argc, char *argv[] )
{
  int c; char *p; FILE *f;
  R lb, ub; /* lower and upper bound of box constraints */
  int option_lubox = 0;
  int iu, io; /* length of proj.transform related vectors in file */
  int out_1T = 0;
  struct cfit_context ctx = 
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0};
  struct cfit_context *C = &ctx;

  /* defaults */
  C->verbose = 0;
  C->maxiter = 60;
  C->frobratio = 1e-6;
  C->QRtolerance = 1e-6;
  C->gamma = -1; /* no projective transform */
  C->chopquantileL = 0; C->chopquantileH = 100;
  
  for(;;)
    {
    int option_index = 0;
    static struct option long_options[] =
        {
          {"help",0,0,'h'},
          {"verbose",0,0,'v'},
          {"Minit",1,0,'i'},
          {"dump",1,0,'d'},
          {"outM",1,0,'M'},
          {"outX",1,0,'X'},
          {"alpha",1,0,'a'},
          {"maxiter",1,0,'T'},
          {"frobratio",1,0,'f'},
          {"constraints",1,0,'C'},
          {"lubox",1,0,'b'},
          {"gamma",1,0,'g'},
          {"viewdir",1,0,'u'},
          {"origin",1,0,'o'},
          {"denom_weight",0,0,'D'},
          {"robust",1,0,'R'},
          {"chop",1,0,'q'},
        
          {"outPM",1,0,1},
          {"outPY",1,0,2},
          {"u_1T",0,0,3},
          {"QRtolerance",1,0,4},
          
          {0,0,0,0}
        };
    c = getopt_long ( argc, argv, "hvk:Y:i:d:M:X:a:T:f:Q:C:x:g:u:o:D:R:q:",
          long_options, &option_index);
    if ( c == -1 ) break;
    switch ( c )
      {
      case 'v':
        C->verbose = 1; break;
      case 'd':
        C->dump_mode = strtol( optarg, &p, 10 );
        if( !p ) C->dump_mode = 0;
        break;
      case 'i':
        C->inM = fopen( optarg,"r" );
        if( ! C->inM ) { perror("Minit"); exit(1); }
        break;
      case 'M':
        if( strcmp(optarg,"stdout") == 0 ) /* allow dumping to stdout */
          {
          if( C->outPM == stdout )
            fprintf(stderr,"Warning: both M and PM ouput to stdout\n");
          C->outM = stdout;
          }
        else
          {
          C->outM = fopen( optarg,"w" );
          if( ! C->outM ) { perror("outM"); exit(1); }
          }
        break;
      case 'X':
        C->outX = fopen( optarg,"w" );
        if( ! C->outX ) { perror("outX"); exit(1); }
        break;
      case 'a':
        {
        int i = 0;
        C->alpha = new_v ( strlen(optarg) );
        p = strtok ( optarg, ", " );
        if( !p )
          { fprintf(stderr,"error in --alpha arguments\n"); exit(1); }
       do 
          {
          if( 1 != sscanf(p,"%"Rf,& C->alpha[i+1] ) ) break;
          i++;
          if( C->alpha[i] < 0 ) C->alpha[i] = 0;
          if( C->alpha[i] > .5 ) C->alpha[i] = .5;
          }
        while ( NULL != (p = strtok ( NULL, ", " )) );
        C->nalpha = i;
        }
        break;
      case 'T':
        C->maxiter = strtol(optarg,&p,10);
        if( p == NULL )
          { fprintf(stderr,"error specifying maxiter\n"); exit(1); }
        break;
      case 'f':
        C->frobratio = strtod(optarg,&p);
        if( p == NULL )
          { fprintf(stderr,"error specifying frobratio\n"); exit(1); }
        break;
      case 'C':
        C->constraints = strtol(optarg,&p,10);
        if( !p || C->constraints < 0 || C->constraints > 3 )
          { fprintf(stderr,"error specifying constraints\n"); exit(1); }
        switch ( C->constraints )
          {
          case 0: C->cf_function = NULL; break;
          case 1: C->cf_function = cf_nn; break;
          case 2: C->cf_function = cf_box; break;
          case 3: C->cf_function = cf_prob; break;
          }
        break;
      case 'x':
        if( C->constraints == 0 )
          { C->constraints = BOX; C->cf_function = cf_box; }
        else if ( C->constraints != BOX )
          break;
        if( 2 != sscanf ( optarg, "%"Rf",%"Rf, &lb, &ub ))
          { fprintf(stderr,"error specifying lubox\n"); exit(1); }
        option_lubox = 1;
        break;

      case 'g':
        C->gamma = strtod( optarg, &p );
        if(!p) { fprintf(stderr,"error specifying gamma\n"); exit(1); }
        break;
      case 'u':
        f = fopen( optarg, "r" );
        C->u = in_newv ( f, &iu );
        fclose(f);
        break;
      case 'o':
        f = fopen( optarg, "r" );
        C->o = in_newv ( f, &io );
        fclose(f);
        break;
      case 'D':
        C->denominator_as_weight = 1;
        break;
      case 'R':
        C->robustcoeff = strtod ( optarg, &p );
        if( !p ) { fprintf(stderr,"robustcoeff?\n"); exit(1); }
        break;
      case 'q':
        C->chopquantileL = strtod ( optarg, &p );
        if( !p ) { fprintf(stderr,"chop?\n"); exit(1); }
        break;
      case 'Q':
        C->chopquantileH = strtod ( optarg, &p );
        if( !p ) { fprintf(stderr,"chop?\n"); exit(1); }
        break;

      case 1:
        if( strcmp(optarg,"stdout") == 0 ) /* allow dumping to stdout */
          {
          if( C->outM == stdout )
            fprintf(stderr,"Warning: both M and PM ouput to stdout\n");
          C->outPM = stdout;
          }
        else
          {
          C->outPM = fopen( optarg,"w" );
          if( ! C->outPM ) { perror("outPM"); exit(1); }
          }
        break;
      case 2:
        C->fn_PY = (char*) malloc (strlen(optarg)+1);
        strcpy( C->fn_PY, optarg );
        break;
      case 3:
        out_1T = 1;
        break;
      case 4:
        C->QRtolerance = strtod(optarg,&p);
        fprintf(stderr,"C->QRtolerance = %"Rg"\n", C->QRtolerance);
        if( p == NULL )
          { fprintf(stderr,"error specifying QRtolerance\n"); exit(1); }
        break;

      case 'h':
      default:
        fprintf(stderr,"%d\n", c);
        usage();
      }
    }
  
  if( argc - optind < 1 ) usage();
  C->k = strtol(argv[optind], &p, 10);
  if( p == NULL )
    {
    fprintf(stderr,"please specify k as the first argument.\n");
    exit(1);
    }
  if( argv[optind+1] )
    {
    C->inY = fopen( argv[optind+1], "r" );
    if( !C->inY )
      { perror(argv[2]); exit(1); }
    }
  else
    C->inY = stdin;

  if( C->nalpha == 0 ) /* default alpha */
    {
    C->alpha = new_v(2);
    C->nalpha = 2;
    C->alpha[1] = 0.5; C->alpha[2] = 0.001;
    }
  
  /*
   * Read data file and initialize matrices
   *
   */
  C->Y = in_newA ( C->inY, &C->n, &C->m );
  if( C->Y ) { fclose( C->inY ); C->inY = NULL; }
  else { fprintf(stderr,"error reading data\n"); exit(1); }
  if(C->verbose)
    fprintf(stderr,"data dimensions: n = %d, m = %d, k = %d\n\n",
        C->n, C->m, C->k );
  
  /* read weight file if specified */
  C->weight = new_v ( C->n );
  if( argv[optind+2] )
    {
    f = fopen( argv[3], "r" );
    if(!f) { perror("weight file"); exit(1); }
    in_v ( f, C->n, C->weight );
    fclose(f);
    }
  else
    c_v ( C->n, 1.0/C->n, C->weight );
  
  /* chop the cone: give a point zero weight if u^T y >= threshold
   *
   * The threshold is specified as a quantile (e.g. include 75% of the points)
   */
  if( C->chopquantileL != 0 || C->chopquantileH != 100 )
    {
    int n = C->n, m = C->m, j;
    R *uTy = new_v ( n ), *suTy = new_v (n), thresholdL, thresholdH;
    if( ! C->u )  /* use the same direction vector as proj transform */
      {
      C->u = new_v ( m );
      if ( out_1T ) c_v ( m, 1, C->u );
      else
        {
        ATv_u ( n, m, C->Y, C->weight, C->u );
        cv ( m, 1.0 / vTv( m, C->u ), C->u );
        }
      }
    for ( j = 1; j <= n; j++ )
      uTy[j] = uTv ( m, C->Y[j], C->u );
    u_v ( n, uTy, suTy );
    qsort ( & suTy[1], n, sizeof(R), &cmpR_asc );
    thresholdL = suTy[ (int) rint(C->chopquantileL * n /100.0) ];
    thresholdH = suTy[ (int) rint(C->chopquantileH * n /100.0) ];
    for ( j = 1; j <= n; j++ )
      if( uTy[j] > thresholdH || uTy[j] < thresholdL )
        C->weight[j] = 0;
    vF1 ( n, C->weight );
    v_del ( uTy );
    v_del ( suTy );
    }
 
  /* projective transform */
  if( C->gamma >= 0 ) /* if gamma negative, ignore */
    {
    int n = C->n, m = C->m, j;
    if( (C->u && iu != m) || (C->o && io != m) )
      {
      fprintf(stderr,"dimensions of vector u and o doesn't match data\n");
      exit(1);
      }
    if( ! C->o )
      {
      C->o = new_v ( m );
      c_v ( m, 0, C->o );
      }
    else
      Am1vT ( n, m, C->Y, C->o );
      
    if( ! C->u )
      {
      C->u = new_v ( m );
      if( out_1T ) c_v ( m, 1, C->u );
      else
        {
        ATv_u ( n, m, C->Y, C->weight, C->u );
        cv ( m, 1.0 / vTv( m, C->u ), C->u );
        }
      }

    for ( j = 1; j <= n; j++ )
      {
      R denom = vFug ( m, C->Y[j], C->u, C->gamma );
      if( C->denominator_as_weight ) C->weight[j] = denom;
      }
    if( C->denominator_as_weight )
      vF1 ( C->n, C->weight );     /* normalize weights */
    }
  
  if( C->fn_PY ) /* output even if gamma < 0, this simply copies Y */
    {
    f = fopen ( C->fn_PY, "w" );
    A_out ( C->n, C->m, C->Y, f );
    fclose(f);
    }


  C->M = new_A ( C->k, C->m );
  C->X = new_A ( C->n, C->k );
  C->Rho = new_A ( C->n, C->m );
  C->normq = new_v ( C->k );
  C->permut = new_P ( C->k );
  C->Mnext = new_A ( C->k, C->m );
  C->Yi = new_A ( C->n, C->m );
  C->Mi = new_A ( C->k, C->m );
  C->W = new_A ( C->n, C->k );
  C->Theta = new_A ( C->k, C->k );
  C->sumWPsi = new_A ( C->k, C->k );
  C->sumWj = new_v ( C->k );

  /* if box constraints is specified, fill in the parameters */
  if( C->constraints == BOX )
    {
    R *p; int m = C->m;
    p = new_v ( 2 * m );
    C->cf_parameters = (void*)p;
    if ( option_lubox )
      { c_v ( m, lb, p ); c_v ( m, ub, p + m ); }
    else
      f = fopen ( "lubox", "r" );
      if( !f )
        {
        fprintf(stderr,"please specify parameters for box constraints\n");
        exit(1);
        }
      in_v ( f, 2 * m, p );
      fclose(f);
    }

  /* setup the initial simplex */
  if( C->inM )
    { in_A ( C->inM, C->k, C->m, C->M ); fclose(C->inM); C->inM = 0; }
  else
    initial_simplex ( C );

  cfit ( C ); /* the damn thing */

  if( C->dump_mode >= FINAL_ONLY )
    cfit_dump ( C );

  exit(0);
}

