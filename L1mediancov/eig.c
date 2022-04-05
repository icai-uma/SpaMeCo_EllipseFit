#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "eig.h"

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg * sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau)

double pythag(double a, double b)
     /* Computes sqrt(a^2 + b^2) without destructive underflow or
	overflow */
{
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);
     
  if (absa > absb) 
    return (absa * sqrt(1.0 + SQR(absb / absa)));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa/absb)));
}


void jacobi (double **a, int n, double *d, double **v, int *nrot)

     /* Computes all eigenvalues and eigenvectors of a real symmetric
	matrix a[1..n][1..n]. On output, elements of a above thep
	diagonal are destroyed.  d[1..n] returns the eigenvalues of
	a.  v[1..n][1..n] is a matrix whose columns contain, on output,
	the normalized eigenvectors of a.  nrot returns the number of
	Jacobi rotations that were required. */
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

  b = (double *) mxCalloc (n, sizeof(double));
  if (b == NULL) {
    mexErrMsgTxt ("mxCalloc b in jacobi()");
    return;
  }
  b--;

  z = (double *) mxCalloc (n, sizeof(double));
  if (z == NULL) {
    mexErrMsgTxt ("mxCalloc z in jacobi()");
    return;
  }
  z--;
 
  /* Initialize to the identity matrix */

  for (ip = 1; ip <= n; ip++) {
    for (iq = 1; iq <= n; iq++) 
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }


  /* Initialize b and d to the diagonal of a.  This vector will
     accumulate terms of the form ta_{pq} as in equation (11.1.14). */

  for (ip = 1; ip <= n; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }
  *nrot = 0;

  for (i=1;i<=50;i++) {
    sm = 0.0;

    /* Sum off-diagonal elements */

    for (ip = 1; ip <= n-1; ip++) {
      for (iq = ip + 1; iq <= n; iq++) 
	sm += fabs(a[ip][iq]);
    }

    /* The normal return, which relies on quadratic convergence to
       machine underflow */

    if (sm == 0.0) {
      mxFree(++z);
      mxFree(++b);
      return;
    }

    if (i < 4)
      tresh = 0.2 * sm / (n*n);	/* on the first three swaps */
    else
      tresh = 0.0;		/* thereafter */


    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1 ; iq<=n; iq++) {
	g = 100.0 * fabs(a[ip][iq]);

	/* After four sweeps, skip the rotation if the off-diagonal
	   element is small. */
	

	if (i > 4 && (double) (fabs(d[ip]) + g) == (double)
	    fabs(d[ip]) && (double) (fabs(d[iq]) + g) == (double)
	    fabs(d[iq]))
	  a[ip][iq] = 0.0;

	    /* Page 2 */

	else if (fabs(a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if ((double) (fabs(h) + g) == (double) fabs(h))
	    t = (a[ip][iq]) / h;	/* t = 1/(2*theta) */
	  else {
	    theta = 0.5 * h / (a[ip][iq]); /* equation 11.1.10 */
	    t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
	    if (theta < 0.0)
	      t = -t;
	  }

	  c = 1.0 / sqrt(1+t*t);
	  s = t * c;
	  tau = s / (1.0 + c);
	  h = t * a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;

	  for (j = 1; j <= ip - 1; j++) { /* Case of rotations
					     1 <= j < p */
	    ROTATE (a, j, ip, j, iq);
	  }
	  for (j = ip + 1; j <= iq - 1; j++) { /* Case of rotations
						  p < j < q */
	    ROTATE (a, ip, j, j, iq);
	  }
	  for (j = iq + 1; j <= n; j++) { /* Case of ratations
					   q < j <= n */
	    ROTATE (a, ip, j, iq, j);
	  }
	  for (j = 1; j <= n; j++) {
	    ROTATE (v, j, ip, j, iq);
	  }
	  ++(*nrot);
	}
      }
    }

    /* Update d with the sum of ta_{pq} and reinitialize z */

    for (ip = 1; ip <= n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  mexErrMsgTxt ("Too many iterations in routine jacobi\n");
}

void eigsrt (double *d, double **v, int n)
     /* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n]
	as output from jacobi (section 11.1) or tqli (section 11.3),
	this routine sorts the eigenvalues into decending order, and
	rearranges the columns of v corespondingly.  The method is
	straight insertion. */
{
  int k, j, i;
  double p;

  for (i = 1; i < n; i++) {
    p = d[k=i];
    for (j = i + 1; j <= n; j++)
      if (fabs(d[j]) >= fabs(p)) 
	p = d[k=j];

    if (k != i) {
      d[k] = d[i];
      d[i] = p;

      for (j = 1; j <= n; j++) {
	p = v[j][i];
	v[j][i] = v[j][k];
	v[j][k] = p;
      }
    }

  }
}

void tred2(double **a, int n, double *d, double *e)
     /* Householder reduction of a real, symmetric matrix
	a[1..n][1..n].  On output, a is replaced by the orthgonal
	matrix Q effecting the transformation.  d[1..n] returns the
	diagonal elements of the tridiagonal matrix, and e[1..n] the
	off-diagonal elements, with e[1] = 0.  Several statements, as
	noted in commensts, can be omitted if only eigenvalues are to
	be found, in which case a contains no useful information on
	output.  Otherwise they are to be included. */
{
  int l, k, j, i;
  double scale, hh, h, g, f;

  for (i = n; i>= 2; i--) {
    l = i - 1;
    h = scale = 0.0;
    if (l > 1) {
      for (k = 1; k <= l; k++) 
	scale += fabs(a[i][k]);
      if (scale == 0.0)		/* skip transformation */
	e[i] = a[i][l];
      else {
	for (k = 1; k <= l; k++) {
	  a[i][k] /= scale;	/* use scaled a's for transformation*/
	  h += a[i][k] * a[i][k]; /* form sigma in h */
	}
	f = a[i][l];
	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i] = scale * g;
	h -= f * g;		/* Now h is equation (11.2.4) */
	a[i][l] = f-g;		/* Store u in the ith row of a. */
	f = 0.0;
	for (j = 1; j <= l; j++) {

	  /* Next statement can be omitted if eigenvectors not wanted
	   */

	  a[j][i] = a[i][j] / h; /* Store u/H in ith column of a. */
	  g = 0.0;		/* Form an element of Au in g. */
	  for (k = 1; k <= j; k++) 
	    g += a[j][k] * a[i][k];
	  for (k = j+1; k <= l; k++)
	    g += a[k][j] * a[i][k];
	  e[j] = g/h;		/* Form element of p in temporarily
				   unused element of e */

	  /* Page 2 */

	  f += e[j] * a[i][j];
	}

	hh = f / (h + h);	/* Form K, equation (11.2.11). */
	for (j = 1; j <= l; j++) { /* Form q and store in e
				      overwriting p */
	  f = a[i][j];		/* Note that e[l] = e[i-1] survives */
	  e[j] = g = e[j] - hh * f;
	  for (k = 1; k <= j; k++) /* Reduce a, equation (11.2.13) */
	    a[j][k] -= (f * e[k] + g * a[i][k]);
	}
      }
    } else
      e[i] = a[i][l];
    d[i] = h;
  }

  /* Next statement can be omitted if eigenvectors not wanted */

  d[1] = 0.0;
  e[1] = 0.0;
  
  /* Contents of this loop can be omitted if eigenvectors not wanted 
   except for statement d[i] = a[i][i]; */

  for (i = 1; i <= n; i++) {	/* Begin accumulation of
				   transformation matrices */
    l = i - 1;
    if (d[i]) {			/* This block skipped when i = 1 */
      for (j = 1; j <= l; j++) {
	g = 0.0;
	for (k = 1; k <= l; k++) /* Use u and u/H stored in a to form
				    PQ */
	  g += a[i][k] * a[k][j];
	for (k = 1; k <= l ; k++)
	  a[k][j] -= g * a[k][i];
      }
    } 
    d[i] = a[i][i];		/* This statement remains */
    a[i][i] = 1.0;		/* Reset row and column of a to
				   identity matrix for next iteration
				*/
    for (j = 1; j <= l; j++)
      a[j][i] = a[i][j] = 0.0;
  }
}

void
tqli(double *d, double *e, int n, double **z)
     /* QL algorithm with implicit shitfs, to determine the
	eigenvalues and eigenvectors of a real, symmetric, tridiagonal
	matrix, or of a real, symmetric matrix previously reduced by
	tred2 (section 11.2).  On input, d[1..n] contains the diagonal
	elements of the tridiagonal matrix.  On output, it returns the
	eigenvalues.  The vector e[1..n] inputs the subdiagonal
	elements of the tridiagonal matrix, with e[1] arbitrary.  On
	output, e is destroyed.  When finding only the eigenvalues,
	several lines may be omitted, as noted in the comments.  If
	the eigenvectors of a tridiagonal matrix are desired, the
	matrix z[1..n][1,,n] is input as the identity matrix.  If the
	eigenvectors of a matrix that has been reduced by tred2 are
	required, then z is input as the matrix output by tred2.  In
	either case, the kth column of z returns the normalized
	eigenvector corresponding to d[k]. */
{
  double pythag (double a, double b);
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  /* Convenient to renumber the elements of e */

  for (i = 2; i <= n; i++) 
    e[i-1] = e[i];
  e[n] = 0.0;
  for (l = 1; l <= n; l++) {
    iter = 0;

    do {

      /* Look for a single small subdiagonal element to split the
	 matrix */

      for (m = l; m <= n - 1; m++) {
	dd = fabs(d[m]) + fabs(d[m+1]);
	if ((double) fabs(e[m] + dd) == dd) 
	  break;
      }

      if (m != l) {
	if (iter++ == (4*n*n+30))
	  mexErrMsgTxt ( "Too many iterations in tqli\n" );
	g = (d[l+1] - d[l]) / (2.0 * e[l]); /* Form shift */
	r = pythag(g, 1.0);
	g = d[m] - d[l] + e[l] / (g + SIGN(r,g)); /* this is d_{m} -
						     k_{s} */
	s = c = 1.0;
	p = 0.0;

	/* Page 2 */

	/* A plane rotation as in the original QL, followed by Givens
	   rotations to restore tridiagonal form. */

	for (i = m-1; i >= l; i--) {

	  f = s * e[i];
	  b = c * e[i];
	  e[i+1] = (r = pythag(f,g));

	  /* recover from underflow */

	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }

	  s = f/r;
	  c = g/r;
	  g = d[i+1] - p;
	  r = (d[i] - g) * s + 2.0 * c * b;
	  d[i+1] = g + (p = s * r);
	  g = c * r - b;

	  /* Next loop can be omitted if eigenvectors not wanted */

	  /* Form eigenvectors */

	  for (k = 1; k <= n; k++) {
	    f = z[k][i+1];
	    z[k][i+1] = s * z[k][i] + c * f;
	    z[k][i] = c * z[k][i] - s * f;
	  }
	}

	if (r == 0.0 && i >= l)
	  continue;

	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}
	  


