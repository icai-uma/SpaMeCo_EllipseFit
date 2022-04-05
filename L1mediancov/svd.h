#ifndef _SVD_H

#define _SVD_H

/*Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
U*W*V'. The matrix U replaces a on output. The diagonal matrix of singular values W is output
as a vector w[1..n]. The matrix V (not the transpose V' ) is output as v[1..n][1..n].*/

void svdcmp(double **a, int m, int n, double w[], double **v);


#endif


