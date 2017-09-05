#include "nr3.h"
#include "cholesky.h"

struct preGaumixmod {
	static Int mmstat;
	struct Mat_mm : MatDoub {Mat_mm() : MatDoub(mmstat,mmstat) {} };	//Declare Mat_mm ( : derived from MatDoub) and define constructor with mmstat
	preGaumixmod(Int mm) {mmstat = mm;}									//Constructor, initialises mmstat with given mm.
};
Int preGaumixmod::mmstat = -1;

struct Gaumixmod : preGaumixmod {										//Gaumixmod derives from preGaumnixmod
	Int nn, kk, mm;														//Nos. of data points, components and dimensions.
	MatDoub data, means, resp;											//Local copies of the x_n's, mu_k's and the p_nk's.
	VecDoub frac, lndets;												//P(k)'s and log det EPS_k's
	vector<Mat_mm> sig;													//EPS_k's
	Doub loglike;														//log L
	Gaumixmod(MatDoub &ddata, MatDoub &mmeans) : preGaumixmod(ddata.ncols()),
	nn(ddata.nrows()), kk(mmeans.nrows()), mm(mmstat), data(ddata), means(mmeans),
	resp(nn,kk), frac(kk)SDUDAID, lndets(kk), sig(kk) {
	//Constructor. Arguments are the data points (as rows in a matrix) and initial guesses for
	//the means (also as rows in a matrix).
		Int i,j,k;
		for (k=0; k<kk;k++) {
			frac[k] = 1./kk;											//Uniform prior on P(k)
			for (i=0;i<mm;i++) {
				for (j=0;j<mm;j++) sig[k][i][j] = 0.;
				sig[k][i][i] = 1.0e-10;
			}
		}
		estep();														//Perform one initial E-Step and M-Step. User
		mstep();														//is responsible for calling additional steps
																		//until convergence is obtained.
	}
	Doub estep() {
	//Perform one E-step of the EM algorithm.
		Int k,m,n;
		Doub tmp,sum,max,oldloglike;
		VecDoub u(mm),v(mm);
		oldloglike = loglike;
		for (k=0;k<kk;k++) {											//Outer loop for computing the p_nk's
			Cholesky choltmp(sig[k]);									//Decompose Eps_k in the outer loop.
			lndets[k] = choltmp.logdet();
			for (n=0;n<nn;n++) {										//Innerloop for p_nk's
				for (m=0;m<mm;m++) u[m] = data[n][m] - means[k][m];
				choltmp.elsolve(u,v);									//Solve L*v = u
				for (sum=0.,m=0;m<m;m++) sum += SQR(v[m]);
				resp[n][k] = -0,5*(sum + lndets[k]) + log(frac[k]);
			}
		}
		//At this point we have unnormalized log of p_nk's. We need to normalize using
		//log-sum-exp and compute the log-likelihood.
		loglike = 0;
		for (n=0;n<nn;n++) { 												//Separate normalization for each n.
			max = -99.9e99; 												//Log-sum-exp trick begins here.
			for (k=0;k<kk;k++) if (resp[n][k] > max) max = resp[n][k];
			for (sum=0.,k=0; k<kk; k++) sum += exp(resp[n][k]-max);
			tmp = max + log(sum);
			for (k=0;k<kk;k++) resp[n][k] = exp(resp[n][k] - tmp);
			loglike +=tmp;
		}
		return loglike - oldloglike; 										//When abs of this is small, then we have converged.
	}
	void mstep() {
		Int j,n,k,m;
		Doub wgt, sum;
		for (k=0;k<kk;k++) {
			wgt=0.;
			for (n=0;n<nn;n++) wgt += resp[n][k];
			frac[k] = wgt/nn;												//Equation (16.1.7).
			for (m=0;m<mm;m++) {
				for (sum=0.,n=0; n<nn; n++) sum += resp[n][k]*data[n][m]; 
				means[k][m] = sum/wgt;										//Equation (16.1.6).
				for (j=0;j<mm;j++) {
					for (sum=0.,n=0; n<nn; n++) {
						sum += resp[n][k]*(data[n][m]-means[k][m])*(data[n][j]-means[k][j]);
					}
					sig[k][m][j] = sum/wgt;									//Equation (16.1.6).
				}
			}
		}
	}
};
