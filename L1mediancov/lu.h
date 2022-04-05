#ifndef _LU_H

#define _LU_H

/* Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation effected by the partial
pivoting; d is output as 1 or -1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.*/
void ludcmp(double **a, int n, int *indx, double *d);

/*Solves the set of n linear equations A*X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, and indx are not modified by this routine
and can be left in place for successive calls with different right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion.*/
void lubksb(double **a, int n, int *indx, double b[]);


/* Hallar la matriz inversa de a[1..n][1..n] mediante la descomposición LU. 
La inversa se devuelve en y[1..n][1..n], mientras que la matriz original a es 
destruida. */
void inv_lu(double **a,double **y,int n);


#endif

