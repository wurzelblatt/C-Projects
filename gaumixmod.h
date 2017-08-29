struct preGaumixmod {
	static Int mmstat;
	struct Mat_mm : MatDoub {Mat_mm() : MatDoub(mmstat,mmstat) {} };	//Declare Mat_mm ( : derived from MatDoub) and define constructor with mmstat
	preGaumixmod(Int mm) {mmstat = mm;}		//Constructor, initialises mmstat with given mm.
};
Int preGaumixmod::mmstat = -1;

struct Gaumixmod : preGaumixmod {		//Gaumixmod derives from preGaumnixmod
	Int nn, kk, mm;						//Nos. of data points, components and dimensions.
	MatDoub data, means, resp;			//Local copies of the x_n's, mu_k's and the p_nk's.
	VecDoub frac, lndets;				//P(k)'s and log det EPS_k's
	vector<Mat_mm> sig;					//EPS_k's
	Doub loglike;						//log L
	Gaumixmod(MatDoub &ddata, MatDoub &mmeans) : preGaumixmod(ddata.ncols()),
	nn(ddata.nrows()), kk(mmeans.nrows()), mm(mmstat), data(ddata), means(mmeans),
	resp(nn,kk), frac(kk), lndets(kk), sig(kk) {
	//Constructor. Arguments are the data points (as rows in a matrix) and initial guesses for
	//the means (also as rows in a matrix).
		Int i,j,k;
		for (k=0; k<kk;k++) {
			frac[k] = 1./kk;				//Uniform prior on P(k)
			for (i=0;i<mm;i++) {
				for (j=0;j<mm;j++) sig[k][i][j] = 0.;
				sig[k][i][i] = 1.0e-10;
			}
		}
		estep();							//Perform one initial E-Step and M-Step. User
		mstep();							//is responsible for calling additional steps
											//until convergence is obtained.
	}
	Doub estep() {
	//Perform one E-step of the EM algorithm.
		Int k,m,n;
		Doub tmp,sum,max,oldloglike;
		VecDoub u(mm),v(mm);
		oldloglike = loglike;
		for (k=0;k<kk;k++) {
			Cholesky choltmp(sig[k]);
			lndets[k] = choltmp.logdet();
			for (n=0;n<nn;n++) {
				for (m=0;m<mm;m++) u[m] = data[n][m] - means[k][m];
				choltmp.esolve(u,v);
				for (sum=0.,m=0;m<m;m++) sum += SQR(v[m]);
				resp[n][k] = -0,5*(sum + lndets[k]) + log(frac[k]);
