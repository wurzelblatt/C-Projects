#include "nr3.h"

struct Cholesky { 																	//cholesky.h
	//Object for Cholesky decomposition of a matrix A, and related functions.
	Int n;
	MatDoub el; 																	//Stores the decomposition.
	Cholesky(MatDoub_I &a) : n(a.nrows()), el(a) {
	//Constructor. Given a positive-definite symmetric matrix a[0..n-1][0..n-1], construct
	//and store its Cholesky decomposition, A = L * L^T .
		Int i,j,k;
		VecDoub tmp;
		Doub sum;
		if (el.ncols() != n) throw("need square matrix");
		for (i=0;i<n;i++) {
			for (j=i;j<n;j++) {
				for (sum=el[i][j],k=i-1;k>=0;k--) sum -= el[i][k]*el[j][k];
				if (i == j) {
					if (sum <= 0.0) 												//A, with rounding errors, is not positive-definite.				
						throw("Cholesky failed");
					el[i][i]=sqrt(sum);
				} else el[j][i]=sum/el[i][i];
			}
		}
		for (i=0;i<n;i++) {
			for (j=0;j<i;j++) {
				el[j][i] = 0.;
			}
		}
	}
	
	void elsolve(VecDoub_I &b, VecDoub_O &y) {
	//Solve L * y = b, where L is the lower triangular matrix in the stored Cholesky decomposition.
	//b[0..n-1] is input as the right-hand side vector.
	//The solution vector is returned in y[0..n-1].
		Int i,j;
		Doub sum;
		if (b.size() != n || y.size() != n) throw("bad lengths");
		for (i=0;i<n;i++) {
			for (sum=b[i],j=0; j<i; j++) sum -= el[i][j]*y[j];
			y[i] = sum/el[i][i];
		}
	}

	Doub logdet() {
	//Return the logarithm of the determinant of A, the matrix whose Cholesky decomposition
	//has been stored.
		Doub sum = 0.;
		for (Int i=0; i<n; i++) sum += log(el[i][i]);
		return 2.*sum;
	}
};
