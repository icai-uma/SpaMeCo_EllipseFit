#ifndef _EIG_H

#define _EIG_H

     /* Computes all eigenvalues and eigenvectors of a real symmetric
	matrix a[1..n][1..n]. On output, elements of a above thep
	diagonal are destroyed.  d[1..n] returns the eigenvalues of
	a.  v[1..n][1..n] is a matrix whose columns contain, on output,
	the normalized eigenvectors of a.  nrot returns the number of
	Jacobi rotations that were required. */

void jacobi (double **a, int n, double *d, double **v, int *nrot);

     /* Householder reduction of a real, symmetric matrix
	a[1..n][1..n].  On output, a is replaced by the orthgonal
	matrix Q effecting the transformation.  d[1..n] returns the
	diagonal elements of the tridiagonal matrix, and e[1..n] the
	off-diagonal elements, with e[1] = 0.  Several statements, as
	noted in commensts, can be omitted if only eigenvalues are to
	be found, in which case a contains no useful information on
	output.  Otherwise they are to be included. */

void tred2(double **a, int n, double *d, double *e);

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

void tqli(double *d, double *e, int n, double **z);

    /* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n]
	as output from jacobi (section 11.1) or tqli (section 11.3),
	this routine sorts the eigenvalues into decending order, and
	rearranges the columns of v corespondingly.  The method is
	straight insertion. */
void eigsrt (double *d, double **v, int n);
     
	
#endif

