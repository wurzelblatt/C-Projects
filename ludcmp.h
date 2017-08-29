#include "nr3.h"


struct LUdcmp
{
	Int n;			
	MatDoub lu;									//Stores the decomposition
	VecInt indx;								//Stores the permutation
	VecDoub vv;									//me
	VecInt P;									//me
	Doub d;										//Used by det
	LUdcmp(MatDoub_I &a);						//Constructor. Argument is the matrix A;
	void solve(VecDoub_I &b, VecDoub_O &x); 	//Solve for a single right-hand side.
	void solve(MatDoub_I &b, MatDoub_O &x);		//Solve for multiple right-hand sides.
	void inverse(MatDoub_O &ainv);				//Calculate matrix inverse A^-1
	void unpermutate(MatDoub_O &aunp);			//unpermutate LUdecomposition; me
	Doub det();									//Return determinant of A
	void mprove(VecDoub_I &b, VecDoub_IO &x);
	MatDoub_I &aref;							//Used only by mprove
};


LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a), aref(a), indx(n)		//Constructor implemented and initialised
{
	const Doub TINY = 1.0e-40;
	Int i, imax,j,k;
	Doub big,temp, tempvv, tempindx;										//me
	//VecDoub vv(n);														//stores implicit scaling of each row.
	vv.resize(n);															//me
	P.resize(n);															//me
	for (i = 0; i < n; i++)													//me
		P[i] = i;														
	d = 1.0;																//No row interchanges yet.
	for (i = 0; i < n; i++)	{												//Loop over rows to get scaling information.
		big = 0.0;
		for (j = 0; j < n; j++) {
			if((temp=abs(lu[i][j])) > big) big=temp;
		}
		if (big == 0.0) throw ("Singular matrix in LUdcmp");				//No nonzero largest element.
		vv[i]=1.0/big;														//Save the scaling.
	}
	for (k = 0; k < n; k++)	{												//Outermost kij loop
		big = 0.0;
		for (i = k; i < n; i++) {
			temp = vv[i] * abs(lu[i][k]);
			//double t1 = vv[i];
			//double t2 = abs(lu[i][k]);
			if (temp > big)	{												//Is the figure of merit for the pivot better than the best so far?
				big = temp;
				imax = i;
			}
		}
		if (k != imax) {													//Do we need to interchange rows?
			for(j = 0; j < n; j++) {										//Yes then do so...
				temp = lu[imax][j];
				lu[imax][j] = lu[k][j];
				lu[k][j] = temp;
			}
			d = -d;															//...and change the parity of d.
			//tempvv = vv[imax];												//me
			vv[imax] = vv[k];												//Also interchange teh scale factor.
			//vv[k] = tempvv;													//me
			tempindx = P[k];												//me
			P[k] = imax;													//me
			P[imax] = tempindx;
		}

		indx[k] = imax;
		
		if (lu[k][k] == 0.0) lu[k][k] = TINY;
		for (i = k+1; i < n; i++) {
			temp = lu[i][k] /= lu[k][k];									//Divide by the pivot element
			
			for (j = k+1; j < n; j++) {										//Innermost loop: reduce remaining submatrix.
				lu[i][j] -= temp*lu[k][j];
			}
		}
	}
}



void LUdcmp::solve(VecDoub_I &b, VecDoub_O &x)
{
	Int i, ii=0, ip, j;
	Doub sum;
	if (b.size() != n || x.size() != n)
		throw("LUdcmp::solve bad sizes");
	for (i=0; i<n;i++) x[i] = b[i];
	for (i=0; i<n;i++) {
		ip=indx[i]; //permutation storage;
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
		else if (sum != 0.0)
			ii = i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum = x[i];
		for (j=i+1;j<n; j++) sum -= lu[i][j]*x[j];
		x[i] = sum/lu[i][i];
	}
}

void LUdcmp::solve(MatDoub_I &b, MatDoub_O &x)
{
	int i,j,m=b.ncols();
	if (b.nrows() != n|| x.nrows() != n || b.ncols() != x.ncols())
		throw("LUdcmp::solve bad sizes");
	VecDoub xx(n);
	for (j=0;j<m;j++) {														//Copy and solve each column in turn.
		for  (i=0;i<n;i++) xx[i] = b[i][j];
		solve(xx,xx);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}


void LUdcmp::inverse(MatDoub_O &ainv)
{
	Int i,j;
	ainv.resize(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv[i][j] = 0.;
		ainv[i][i] = 1.;
	}
	solve(ainv,ainv); 
}

Doub LUdcmp::det()
{
	Doub dd = d;
	for (Int i=0;i<n;i++) dd *= lu[i][i];
	return dd;
}

void LUdcmp::unpermutate(MatDoub_O &aunp)
{
	Doub temp;
	Int k,j;
	aunp.resize(n,n);
	for (k = 0; k < n; k++)	{
		if (k != indx[k]) {
			for(j = 0; j < n; j++) {									
				temp = lu[indx[k]][j];
				lu[indx[k]][j] = lu[k][j];
				lu[k][j] = temp;
			}
		}
	
	}
	
		
}
	

	
			
			
